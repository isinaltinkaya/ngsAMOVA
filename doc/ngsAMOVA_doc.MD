# ngsAMOVA Documentation
## Table of Contents
- [AMOVA Formula](#amova-formula)
- [Metadata File Input](#metadata-file-input)
- [EM Optimization](#em-optimization)
- [AMOVA Analysis](#amova-analysis)
- [Block Bootstrapping](#block-bootstapping)
- [Distance Matrix Output](#distance-matrix-output)
- [Distance Matrix File Input](#distance-matrix-file-input)


## AMOVA Formula 

The formula option `--formula` defines the AMOVA model to be fitted. The formula must be of the form:

```sh
--formula Individual ~ Group/subgroup/subsubgroup
```

The formula can be modified to suit the specific experimental design. 

- The left-hand side of the tilde `~` is the AMOVA model to be fitted with the distance matrix. 
- The right-hand side specifies the hierarchical grouping (population, region, etc). The groups are in the order of decreasing hierarchy and are separated by `/`.
- Specifying `--formula Individual ~ Population` will simply test for population differentiation.
- The group names in the formula must exist among the column names in the metadata file. For details about the metadata file, see [Metadata File Input](#metadata-file-input---metadata).


## Metadata File Input 

The metadata file is specified using the `--metadata <filename>` option. Metadata file is a tab-separated file that contains the information about the samples-to-groups mapping. 

- The metadata file can contain any number of columns. 
- The names of the values listed below the grouping that was used in the left-hand side of the AMOVA formula must match the names of the individuals as they appear in the input data file, e.g. VCF file or distance matrix.
- The metadata file must contain a header row.

For example, the following metadata file contains information about the individuals, regions, populations, and subpopulations:

```
$ cat tests/data/metadata_Individual_Region_Population_Subpopulation.tsv
Individual    Region    Population    Subpopulation
pop1_ind1     reg1      pop1          subpop1
pop1_ind2     reg1      pop1          subpop2
pop1_ind3     reg1      pop1          subpop2
pop2_ind1     reg1      pop2          subpop3
pop2_ind2     reg1      pop2          subpop3
pop2_ind3     reg1      pop2          subpop3
pop3_ind1     reg2      pop3          subpop4
pop3_ind2     reg2      pop3          subpop4
pop3_ind3     reg2      pop3          subpop4
```

Here, if a formula `Individual~Region/Population/Subpopulation` was used, we can read this metadata file as "the Individual `pop1_ind1` belongs to the Subpopulation `subpop1`, which belongs to the Population `pop1` in the Region `reg1`".

## EM Optimization

Using `--verbose 1` will print the reason for termination.

## AMOVA Analysis

- **AMOVA Output (`--print-amova`)**
	-  `--print-amova 1` will print the AMOVA results in CSV format (`output.amova.csv`)
	- `--print-amova 2` will print the AMOVA results in table format (`output.amova_table.txt`)

Example AMOVA analysis:

```sh
./bin/ngsAMOVA \
	--in-vcf ./tests/data/test_s9_d1_1K_2contigs.vcf \
	--bcf-src 1 \      # Use genotype likelihoods (GL field in the VCF file)
	-doMajorMinor 1 \  # Use REF and first ALT as major and minor 
	-doEM 1 \          # Use EM algorithm (necessary to do -doJGTM with GL data i.e. --bcf-src 1)
	-doJGTM 1 \        # Estimate the joint genotypes matrix (necessary for -doDist)
	-doDist 1 \        # Calculate pairwise distances (necessary for -doAMOVA)
	-doAMOVA 1 \       # Perform AMOVA analysis
	--print-amova 3 \  # Print AMOVA results in both CSV and table format
	--minInd 2 \       # Minimum number of individuals with data for a site to be included in analyses
	--maxEmIter 10 \   # Maximum number of EM iterations (termination criteria)
	--metadata ./tests/data/metadata_Individual_Population.tsv \
	--formula 'Individual~Population'
```

The metadata file `metadata_Individual_Population.tsv` is a tab-separated values file containing the individual-to-population mapping. For details, see the section [Metadata file input](#metadata-file-input---metadata).

```
$ cat metadata_Individual_Population.tsv
Individual   Population
pop1_ind1    pop1
pop1_ind2    pop1
pop1_ind3    pop1
pop2_ind1    pop2
pop2_ind2    pop2
pop2_ind3    pop2
pop3_ind1    pop3
pop3_ind2    pop3
pop3_ind3    pop3
```

Here, the metadata file contains two columns: `Individual` and `Population`. Each individual is assigned to a population. For example, `pop1_ind1` belongs to population `pop1`.

The formula `Individual~Population` specifies the AMOVA model. Here, the individuals are nested within the populations. For more details on the formula, see the section [AMOVA formula](#amova-formula---formula).

Running the above command in ngsAMOVA directory will produce the following output files:

```
$ cat output.amova.csv 
df,Total,8
SSD,Total,0.031619158643330654
MSD,Total,0.0039523948304163318
df,Among_Population_within_Total,2
SSD,Among_Population_within_Total,0.010226069566826711
MSD,Among_Population_within_Total,0.0051130347834133557
df,Among_Individual_within_Population,6
SSD,Among_Individual_within_Population,0.021393089076503943
MSD,Among_Individual_within_Population,0.0035655148460839903
Phi,Population_in_Total,0.1263893979336376
Variance_coefficient,c_0_0,3
Variance_coefficient,c_0_1,1
Variance_coefficient,c_1_1,1
Variance_component,Population,0.0005158399791097884
Percentage_variance,Population,12.638939793363759
Variance_component,Individual,0.0035655148460839903
Percentage_variance,Individual,87.361060206636239
```

```
$ cat output.amova_table.txt 
================================================================================
=== AMOVA ======================================================================
Formula: Individual~Population

Source of variation                              d.f.      SSD          MSD
--------------------------------------------------------------------------------
Among Population within Total                    2         0.010226     0.005113     
Among Individual within Population               6         0.021393     0.003566     

Variance coefficients:
c_0_0	3
c_0_1	1
c_1_1	1

Variance components:
Population                    	0.000516  	(12.638940  %)
Individual                    	0.003566  	(87.361060  %)

Phi statistics:
Phi(Population in Total)	0.126389

================================================================================
```

## Block bootstapping

Block bootstrapping test is enabled via `-doBlockBootstrap 1` option. The blocks can be generated in three ways:

- **With `--block-size <int>`**

The genome is divided into blocks of size `--block-size <int>` (default: unset). The blocks are then sampled with replacement to create bootstrap replicates.

- **With `--in-blocks-tab <file>`**

The blocks are read from a tab-separated values file. The file must have 3 columns: contig name, block start position, block end position. The positions are 1-based and inclusive for both start and end.

Example blocks tab file:

```
chr22    16100001    16200000
chr22    16200001    16300000
```

- **With `--in-blocks-bed <file>`**

The blocks are read from a BED file. The file must have 3 columns: contig name, block start position, block end position. Positions are 0-based and inclusive for the start and exclusive for the end, as defined in the BED file standards. 

Example blocks BED file:

```
chr22    16100000    16200001
chr22    16200000    16300001
```

### Printing the blocks

The blocks can be printed to a file using the `--print-blocks 1` option. The blocks file `output.blocks.tsv` has 4 columns: 

- **Contig**: Name of the contig the block belongs to
- **Start**: Start position of the block (1-based inclusive)
- **End**: End position of the block (1-based inclusive)
- **NSites**: Number of sites in the block

Example:

```sh
./bin/ngsAMOVA \
	--in-vcf tests/data/test_s9_d1_1K_2contigs_emptyBlock.bcf \
	--bcf-src 1 \
	-doMajorMinor 1 \
	-doJGTM 1 \
	-doEM 1 \
	--print-dm 3 \
	-doDist 1 \
	--maxEmIter 10 \
	-doBlockBootstrap 1 \
	--print-blocks 1 \
	--minInd 9 \
	--block-size 100000 \
	--nBootstraps 2 \
	--seed 3
```

```
$ cat output.blocks.tsv
Contig   Start	End	       NSites
chr22    16000001	16100000   4
chr22    16100001	16200000   0
chr23    6100001	6200000    4
chr23    6200001	6300000    5
```

 Here, notice that the second block of contig `chr22` has 0 sites. This is because there are no sites left in this block after applying filters (e.g. `--minInd 9`). This block will still be included in the bootstrap sampling.

### Printing the bootstrap samples

The bootstrap samples can be printed to a file using the `--print-bootstrap 1` option. The bootstrap samples file `output.bootstrap_samples.tsv` has 7 columns: 

 - **Rep**: Replicate number
 - **Pos**: Block position in the genome
 - **BlockID**: Block ID (0-based, equal to the index of the block in the blocks file)
 - **BlockContig**: Contig of the block
 - **BlockStart**: Start position of the block (1-based inclusive)
 - **BlockEnd**: End position of the block (1-based inclusive)
 - **NSites**: Number of sites in the block

Example:

```sh
./bin/ngsAMOVA \
	--in-vcf tests/data/test_s9_d1_1K_2contigs_emptyBlock.bcf \
	--bcf-src 1 \
	-doMajorMinor 1 \
	-doJGTM 1 \
	-doEM 1 \
	--print-dm 3 \
	-doDist 1 \
	--maxEmIter 10 \
	-doBlockBootstrap 1 \
	--minInd 9 \
	--block-size 100000 \
	--print-bootstrap 1 \
	--nBootstraps 2 \
	--seed 3
```

```
$ cat output.bootstrap_samples.tsv
Rep  Pos  BlockID  BlockContig  BlockStart   BlockEnd     NSites
0    0    3        chr23        6200001      6300000      5
0    1    3        chr23        6200001      6300000      5
0    2    1        chr22        16100001     16200000     0
0    3    1        chr22        16100001     16200000     0
1    0    2        chr23        6100001      6200000      4
1    1    1        chr22        16100001     16200000     0
1    2    0        chr22        16000001     16100000     4
1    3    0        chr22        16000001     16100000     4
```

Here, the bootstrap samples are generated by sampling with replacement from the blocks. The number of bootstrap samples is equal to the number of blocks. Each row in the file represents a block in a bootstrap replicate. The `Pos` column represents the position of the block in the replicate's synthetic genome (0-based). 

The genome arrangement in this example, where we have 4 blocks in total for 2 contigs, is:

```
[block]  [block]  [block]  [block]
|_Pos 0  |_Pos 1  |_Pos 2  |_Pos 3
```

The first bootstrap replicate is arranged as follows:

```
[BlockId:3] [BlockId:3] [BlockId:1] [BlockId:1]
```

The second bootstrap replicate is arranged as follows:

```
[BlockId:2] [BlockId:1] [BlockId:0] [BlockId:0]
```

This arrangement corresponds to the following genome for the first bootstrap replicate:

```
[chr23:6200001-6300000] [chr23:6200001-6300000] [chr22:16100001-16200000] [chr22:16100001-16200000]
```

and for the second bootstrap replicate:

```
[chr23:6100001-6200000] [chr22:16100001-16200000] [chr22:16000001-16100000] [chr22:16000001-16100000]
```

You can think of this as cutting these pieces of the genome and pasting them together to create a new synthetic genome to perform the bootstrap analyses on. We perform this the same way for each individual to preserve the correlation structure of the data.

## Distance Matrix Output

Example analysis:

```sh
./bin/ngsAMOVA \
	--in-vcf ./tests/data/test_s9_d1_1K_2contigs.vcf \
	--bcf-src 1 \       # Use genotype likelihoods (GL field in the VCF file)
	-doMajorMinor 1 \   # Use REF and first ALT as major and minor
	-doEM 1 \           # Use EM algorithm (necessary to do -doJGTM with GL data i.e. --bcf-src 1)
	-doJGTM 1 \         # Estimate the joint genotypes matrix (necessary for -doDist)
	-doDist 1 \		    # Calculate pairwise distances 
	--print-dm 3 \      # Print distance matrix in both condensed and verbose format
	--maxEmIter 10      # Maximum number of EM iterations (termination criteria)
```


- **Distance Matrix Output (`--print-dm`)**
	- `--print-dm 1` will print the distance matrix in condensed format (`output.distance_matrix.txt`)
	- `--print-dm 2` will print the distance matrix in verbose format (`output.distance_matrix.verbose.txt`)

- **Distance Matrix in Condensed Format (`--print-dm 1`)**

`--print-dm 1` will print the distance matrix in condensed format: (`output.distance_matrix.txt`)


```
$ cat output.distance_matrix.txt 
1
LTED
0
0
9
36
pop1_ind1
pop1_ind2
pop1_ind3
pop2_ind1
pop2_ind2
pop2_ind3
pop3_ind1
pop3_ind2
pop3_ind3
0.099941876464224796
0.095928932519746588
0.07466008715140468
0.06988508306555366
0.086496181719053961
0.075635796461499746
0.084242255586261458
0.090241760274019894
0.083620725852630529
0.057715408482645086
0.091270317689165714
0.065973961092172989
0.099600316460724259
0.072870937444969391
0.065332554796685954
0.083582025378800842
0.080412582081914638
0.083409695099543088
0.059369524410479324
0.088809269953640038
0.063974727955852406
0.11055160045412826
0.10501214425273189
0.10880269760987354
0.084970426798190987
0.093733459221089208
0.061528175018082616
0.088630642729381981
0.099248349010176609
0.091536067699058407
0.11064248101093817
0.072865827284072787
0.082177413817313902
0.066113281888749556
0.078121591759825998
0.079924397515775986
```

Line number | Description
------------|------------
1           | Number of distance matrices
2           | Type of the distance matrix (2: LTED, 3: LTID, 4: UTED, 5: UTID, 7: FULL)
3 		    | Transformation applied to the distances (0: None, 1: Square)
4           | Distance metric used to calculate the distances (0: Dij, 1: Sij, 2: Fij, 3: IBS0, 4: IBS1, 5: IBS2, 6: Kinship, 7: R0, 8: R1, 9: Dxy)
5           | Number of items in the distance matrix
6           | Number of distances in the distance matrix
[7-X)       | Names of the items in the distance matrix (X = 6 + line 5)
[X-Y)       | Pairwise distances in the distance matrix (Y = X + line 6)



Distance matrix types:

Type | Description
-----|------------
LTED | Lower triangular, excluding the diagonal
LTID | Lower triangular, including the diagonal
UTED | Upper triangular, excluding the diagonal
UTID | Upper triangular, including the diagonal
FULL | Full matrix

Using R, this distance matrix can be read using the following code:

```R
d <- read.csv("output.distance_matrix.txt",header=FALSE)
nInd <- as.numeric(d[5,])                       # Number of individuals
nIndPairs <- as.numeric(d[6,])                  # Number of pairwise distances
indNames <- d[(7):(7+nInd-1),] 				    # Names of the individuals
distances <- d[(7+nInd):(7+nInd+nIndPairs-1),]  # Pairwise distances
distances[distances=="nan"]<-NA 			    # Replace "nan" with NA, if any
```

- **Distance Matrix in Verbose Format (`--print-dm 2`)**

`--print-dm 2` will print the distance matrix in verbose format: (`output.distance_matrix.verbose.txt`)

```
$ cat output.distance_matrix.verbose.txt
## This file was produced by: ngsAMOVA [version: v0.6] [build: Feb 24 2025 15:27:00] [htslib: 1.15.1-33-ge9781647]
## Command: ngsAMOVA --in-vcf /home/isin/Projects/AMOVA/recon/tests/data/data8_1.vcf --bcf-src 1 -doMajorMinor 1 -doBlockBootstrap 1 -doJGTM 1 -doEM 1 --print-dm 1 --minInd 2 -doDist 1 --maxEmIter 10 --rm-invar-sites 1 --seed 2 --block-size 50000 --print-blocks 1 --nbootstraps 1 --allow-mispairs 1 --min-npairs 0 --print-dm 3
## Number of distance matrices: 2
## -> Following information is shared for all distance matrices in this file:
# Type of the distance matrix: Lower Triangular Matrix (Excluding Diagonal)
# Transformation applied to the distances: None (--dm-transform 0)
# Method used to calculate the distances: Dij (--dm-method 0)
# Number of items in the distance matrix: 40
# Number of distances in the distance matrix: 780
# NAME, Names of the items in the distance matrix:
# NAME	[2]INDEX	[3]NAME                  
NAME  	0       	popA_ind1                
NAME  	1       	popA_ind2                
NAME  	2       	popA_ind3                
NAME  	3       	popA_ind4                
NAME  	4       	popA_ind5                
NAME  	5       	popA_ind6                
NAME  	6       	popA_ind7                
NAME  	7       	popA_ind8                
NAME  	8       	popA_ind9                
NAME  	9       	popA_ind10               
NAME  	10      	popB_ind1                
NAME  	11      	popB_ind2                
NAME  	12      	popB_ind3                
NAME  	13      	popB_ind4                
NAME  	14      	popB_ind5                
NAME  	15      	popB_ind6                
NAME  	16      	popB_ind7                
NAME  	17      	popB_ind8                
NAME  	18      	popB_ind9                
NAME  	19      	popB_ind10               
NAME  	20      	popC_ind1                
NAME  	21      	popC_ind2                
NAME  	22      	popC_ind3                
NAME  	23      	popC_ind4                
NAME  	24      	popC_ind5                
NAME  	25      	popC_ind6                
NAME  	26      	popC_ind7                
NAME  	27      	popC_ind8                
NAME  	28      	popC_ind9                
NAME  	29      	popC_ind10               
NAME  	30      	popD_ind1                
NAME  	31      	popD_ind2                
NAME  	32      	popD_ind3                
NAME  	33      	popD_ind4                
NAME  	34      	popD_ind5                
NAME  	35      	popD_ind6                
NAME  	36      	popD_ind7                
NAME  	37      	popD_ind8                
NAME  	38      	popD_ind9                
NAME  	39      	popD_ind10               
## -> Following information is specific to distance matrix with MATRIX_INDEX=0 in this file:
# DIST(MATRIX_INDEX), Pairwise distances in the distance matrix:
# DIST(0)        	[4]PAIR_INDEX	[5]FIRST_ITEM_INDEX	[6]SECOND_ITEM_INDEX	[7]FIRST_ITEM_NAME       	[8]SECOND_ITEM_NAME      	[9]DISTANCE              
DIST(0)          	0            	1                  	0                  	popA_ind2                	popA_ind1                	0.11467617878074979      
DIST(0)          	1            	2                  	0                  	popA_ind3                	popA_ind1                	0.093092202947423874     
DIST(0)          	2            	2                  	1                  	popA_ind3                	popA_ind2                	0.043429845634977288     
DIST(0)          	3            	3                  	0                  	popA_ind4                	popA_ind1                	0.09906854300806188      
DIST(0)          	4            	3                  	1                  	popA_ind4                	popA_ind2                	0.14476271136798249      
DIST(0)          	5            	3                  	2                  	popA_ind4                	popA_ind3                	0.09250928862574137      
DIST(0)          	6            	4                  	0                  	popA_ind5                	popA_ind1                	0.11548749475140634      
DIST(0)          	7            	4                  	1                  	popA_ind5                	popA_ind2                	0.10235826812419881      
DIST(0)          	8            	4                  	2                  	popA_ind5                	popA_ind3                	0.095219871094828198     
DIST(0)          	9            	4                  	3                  	popA_ind5                	popA_ind4                	0.10924768687545443      
DIST(0)          	10           	5                  	0                  	popA_ind6                	popA_ind1                	0.10717339691795179      
DIST(0)          	11           	5                  	1                  	popA_ind6                	popA_ind2                	0.081245857649091494     
DIST(0)          	12           	5                  	2                  	popA_ind6                	popA_ind3                	0.10665280767883598      
DIST(0)          	13           	5                  	3                  	popA_ind6                	popA_ind4                	0.0992253751337005       
DIST(0)          	14           	5                  	4                  	popA_ind6                	popA_ind5                	0.084457781887820244     
DIST(0)          	15           	6                  	0                  	popA_ind7                	popA_ind1                	0.10530027483526966      
DIST(0)          	16           	6                  	1                  	popA_ind7                	popA_ind2                	0.10952915292971824      
DIST(0)          	17           	6                  	2                  	popA_ind7                	popA_ind3                	0.056029316909226953     
DIST(0)          	18           	6                  	3                  	popA_ind7                	popA_ind4                	0.11165936466273685      
DIST(0)          	19           	6                  	4                  	popA_ind7                	popA_ind5                	0.11752068019191365      
DIST(0)          	20           	6                  	5                  	popA_ind7                	popA_ind6                	0.099017757040388107     
DIST(0)          	21           	7                  	0                  	popA_ind8                	popA_ind1                	0.097166560252935227     
DIST(0)          	22           	7                  	1                  	popA_ind8                	popA_ind2                	0.087725183305922882     
DIST(0)          	23           	7                  	2                  	popA_ind8                	popA_ind3                	0.035463864283368375     
DIST(0)          	24           	7                  	3                  	popA_ind8                	popA_ind4                	0.10971550268283765      
DIST(0)          	25           	7                  	4                  	popA_ind8                	popA_ind5                	0.052211145413247576     
DIST(0)          	26           	7                  	5                  	popA_ind8                	popA_ind6                	0.080090237697028069     
DIST(0)          	27           	7                  	6                  	popA_ind8                	popA_ind7                	0.081349577339648454     
DIST(0)          	28           	8                  	0                  	popA_ind9                	popA_ind1                	0.090441910531874212     
DIST(0)          	29           	8                  	1                  	popA_ind9                	popA_ind2                	0.095867900816094878     
DIST(0)          	30           	8                  	2                  	popA_ind9                	popA_ind3                	0.058549180362201261     
DIST(0)          	31           	8                  	3                  	popA_ind9                	popA_ind4                	0.1026648719420762       
DIST(0)          	32           	8                  	4                  	popA_ind9                	popA_ind5                	0.069954678492391426     
DIST(0)          	33           	8                  	5                  	popA_ind9                	popA_ind6                	0.092754438818220819     
DIST(0)          	34           	8                  	6                  	popA_ind9                	popA_ind7                	0.067255283583122355     
DIST(0)          	35           	8                  	7                  	popA_ind9                	popA_ind8                	0.081902410789416374     
DIST(0)          	36           	9                  	0                  	popA_ind10               	popA_ind1                	0.12644992674147426      
DIST(0)          	37           	9                  	1                  	popA_ind10               	popA_ind2                	0.096636012476263314     
DIST(0)          	38           	9                  	2                  	popA_ind10               	popA_ind3                	0.10716580364420708      
DIST(0)          	39           	9                  	3                  	popA_ind10               	popA_ind4                	0.10184357392425492      
DIST(0)          	40           	9                  	4                  	popA_ind10               	popA_ind5                	0.10530742687306804      
DIST(0)          	41           	9                  	5                  	popA_ind10               	popA_ind6                	0.079663031092544523     
DIST(0)          	42           	9                  	6                  	popA_ind10               	popA_ind7                	0.13024246593697633      
DIST(0)          	43           	9                  	7                  	popA_ind10               	popA_ind8                	0.08651233086631463      
DIST(0)          	44           	9                  	8                  	popA_ind10               	popA_ind9                	0.091715136836088701     
DIST(0)          	45           	10                 	0                  	popB_ind1                	popA_ind1                	0.2336180839991483       
DIST(0)          	46           	10                 	1                  	popB_ind1                	popA_ind2                	0.1895475907764064       
DIST(0)          	47           	10                 	2                  	popB_ind1                	popA_ind3                	0.2028394442494629       
DIST(0)          	48           	10                 	3                  	popB_ind1                	popA_ind4                	0.23700842812959011      
DIST(0)          	49           	10                 	4                  	popB_ind1                	popA_ind5                	0.21720537528782002      
DIST(0)          	50           	10                 	5                  	popB_ind1                	popA_ind6                	0.16646532629003258      
DIST(0)          	51           	10                 	6                  	popB_ind1                	popA_ind7                	0.1794056680155674       
DIST(0)          	52           	10                 	7                  	popB_ind1                	popA_ind8                	0.18782183744554481      
DIST(0)          	53           	10                 	8                  	popB_ind1                	popA_ind9                	0.17880978453375057      
DIST(0)          	54           	10                 	9                  	popB_ind1                	popA_ind10               	0.13315769854599507      
DIST(0)          	55           	11                 	0                  	popB_ind2                	popA_ind1                	0.15574930521351982      
DIST(0)          	56           	11                 	1                  	popB_ind2                	popA_ind2                	0.12902711380341364      
DIST(0)          	57           	11                 	2                  	popB_ind2                	popA_ind3                	0.18546996912365243      
DIST(0)          	58           	11                 	3                  	popB_ind2                	popA_ind4                	0.2009515144952744       
DIST(0)          	59           	11                 	4                  	popB_ind2                	popA_ind5                	0.13921199302111537      
DIST(0)          	60           	11                 	5                  	popB_ind2                	popA_ind6                	0.1186383060851682       
DIST(0)          	61           	11                 	6                  	popB_ind2                	popA_ind7                	0.2186548244483244       
DIST(0)          	62           	11                 	7                  	popB_ind2                	popA_ind8                	0.16953188761613003      
DIST(0)          	63           	11                 	8                  	popB_ind2                	popA_ind9                	0.13355896813518409      
DIST(0)          	64           	11                 	9                  	popB_ind2                	popA_ind10               	0.13685061378832791      
DIST(0)          	65           	11                 	10                 	popB_ind2                	popB_ind1                	0.082160885097566222     
DIST(0)          	66           	12                 	0                  	popB_ind3                	popA_ind1                	0.18435734223985237      
DIST(0)          	67           	12                 	1                  	popB_ind3                	popA_ind2                	0.18073172098884421      
DIST(0)          	68           	12                 	2                  	popB_ind3                	popA_ind3                	0.12921894115135346      
DIST(0)          	69           	12                 	3                  	popB_ind3                	popA_ind4                	0.12284745856950652      
DIST(0)          	70           	12                 	4                  	popB_ind3                	popA_ind5                	0.18104606093364706      
DIST(0)          	71           	12                 	5                  	popB_ind3                	popA_ind6                	0.19471725344017071      
DIST(0)          	72           	12                 	6                  	popB_ind3                	popA_ind7                	0.19364555475485268      
DIST(0)          	73           	12                 	7                  	popB_ind3                	popA_ind8                	0.11517994569806497      
DIST(0)          	74           	12                 	8                  	popB_ind3                	popA_ind9                	0.1348604073845501       
DIST(0)          	75           	12                 	9                  	popB_ind3                	popA_ind10               	0.19398220480765638      
DIST(0)          	76           	12                 	10                 	popB_ind3                	popB_ind1                	0.16274535265075507      
DIST(0)          	77           	12                 	11                 	popB_ind3                	popB_ind2                	0.13373797837552712      
DIST(0)          	78           	13                 	0                  	popB_ind4                	popA_ind1                	0.13365602036443922      
DIST(0)          	79           	13                 	1                  	popB_ind4                	popA_ind2                	0.17171434159611629      
DIST(0)          	80           	13                 	2                  	popB_ind4                	popA_ind3                	0.14385416636385348      
DIST(0)          	81           	13                 	3                  	popB_ind4                	popA_ind4                	0.13620392951227267      
DIST(0)          	82           	13                 	4                  	popB_ind4                	popA_ind5                	0.13844854551032743      
DIST(0)          	83           	13                 	5                  	popB_ind4                	popA_ind6                	0.20507900555218073      
DIST(0)          	84           	13                 	6                  	popB_ind4                	popA_ind7                	0.16674800759215752      
DIST(0)          	85           	13                 	7                  	popB_ind4                	popA_ind8                	0.19158365508898345      
DIST(0)          	86           	13                 	8                  	popB_ind4                	popA_ind9                	0.14481723658856765      
DIST(0)          	87           	13                 	9                  	popB_ind4                	popA_ind10               	0.17216895148979988      
DIST(0)          	88           	13                 	10                 	popB_ind4                	popB_ind1                	0.10146661607076186      
DIST(0)          	89           	13                 	11                 	popB_ind4                	popB_ind2                	0.12618862661658062      
DIST(0)          	90           	13                 	12                 	popB_ind4                	popB_ind3                	0.058739516198463057     
DIST(0)          	91           	14                 	0                  	popB_ind5                	popA_ind1                	0.21820072728406747      
DIST(0)          	92           	14                 	1                  	popB_ind5                	popA_ind2                	0.1465067556144577       
DIST(0)          	93           	14                 	2                  	popB_ind5                	popA_ind3                	0.19196571222817477      
DIST(0)          	94           	14                 	3                  	popB_ind5                	popA_ind4                	0.1761893161168778       
DIST(0)          	95           	14                 	4                  	popB_ind5                	popA_ind5                	0.16128517489385935      
DIST(0)          	96           	14                 	5                  	popB_ind5                	popA_ind6                	0.15657767928982227      
DIST(0)          	97           	14                 	6                  	popB_ind5                	popA_ind7                	0.17567679918954104      
DIST(0)          	98           	14                 	7                  	popB_ind5                	popA_ind8                	0.16873546816436635      
DIST(0)          	99           	14                 	8                  	popB_ind5                	popA_ind9                	0.19334619577350509      
DIST(0)          	100          	14                 	9                  	popB_ind5                	popA_ind10               	0.14622659497598675      
DIST(0)          	101          	14                 	10                 	popB_ind5                	popB_ind1                	0.076306207154471978     
DIST(0)          	102          	14                 	11                 	popB_ind5                	popB_ind2                	0.066529877811036531     
DIST(0)          	103          	14                 	12                 	popB_ind5                	popB_ind3                	0.1409409842754705       
DIST(0)          	104          	14                 	13                 	popB_ind5                	popB_ind4                	0.096059387487813133     
DIST(0)          	105          	15                 	0                  	popB_ind6                	popA_ind1                	0.13848934374173905      
DIST(0)          	106          	15                 	1                  	popB_ind6                	popA_ind2                	0.15612309594158841      
DIST(0)          	107          	15                 	2                  	popB_ind6                	popA_ind3                	0.17415526086063793      
DIST(0)          	108          	15                 	3                  	popB_ind6                	popA_ind4                	0.19717244117112806      
DIST(0)          	109          	15                 	4                  	popB_ind6                	popA_ind5                	0.20189920751847495      
DIST(0)          	110          	15                 	5                  	popB_ind6                	popA_ind6                	0.21586937475895537      
DIST(0)          	111          	15                 	6                  	popB_ind6                	popA_ind7                	0.18762449822932767      
DIST(0)          	112          	15                 	7                  	popB_ind6                	popA_ind8                	0.19803684293859386      
DIST(0)          	113          	15                 	8                  	popB_ind6                	popA_ind9                	0.17157542180045349      
DIST(0)          	114          	15                 	9                  	popB_ind6                	popA_ind10               	0.152054208614887        
DIST(0)          	115          	15                 	10                 	popB_ind6                	popB_ind1                	0.1370090321125082       
DIST(0)          	116          	15                 	11                 	popB_ind6                	popB_ind2                	0.087766385293120752     
DIST(0)          	117          	15                 	12                 	popB_ind6                	popB_ind3                	0.14859284705566383      
DIST(0)          	118          	15                 	13                 	popB_ind6                	popB_ind4                	0.13156715752246437      
DIST(0)          	119          	15                 	14                 	popB_ind6                	popB_ind5                	0.10742285964688827      
DIST(0)          	120          	16                 	0                  	popB_ind7                	popA_ind1                	0.17973240429294551      
DIST(0)          	121          	16                 	1                  	popB_ind7                	popA_ind2                	0.14078806644125302      
DIST(0)          	122          	16                 	2                  	popB_ind7                	popA_ind3                	0.091761684712605401     
DIST(0)          	123          	16                 	3                  	popB_ind7                	popA_ind4                	0.14915398801606838      
DIST(0)          	124          	16                 	4                  	popB_ind7                	popA_ind5                	0.16024652644638926      
DIST(0)          	125          	16                 	5                  	popB_ind7                	popA_ind6                	0.16262642621955026      
DIST(0)          	126          	16                 	6                  	popB_ind7                	popA_ind7                	0.12382788967982038      
DIST(0)          	127          	16                 	7                  	popB_ind7                	popA_ind8                	0.12664600605320159      
DIST(0)          	128          	16                 	8                  	popB_ind7                	popA_ind9                	0.11738137106945715      
DIST(0)          	129          	16                 	9                  	popB_ind7                	popA_ind10               	0.15210423771205939      
DIST(0)          	130          	16                 	10                 	popB_ind7                	popB_ind1                	0.082060786862680438     
DIST(0)          	131          	16                 	11                 	popB_ind7                	popB_ind2                	0.13608893161450786      
DIST(0)          	132          	16                 	12                 	popB_ind7                	popB_ind3                	0.080284930300773313     
DIST(0)          	133          	16                 	13                 	popB_ind7                	popB_ind4                	0.069619761483182135     
DIST(0)          	134          	16                 	14                 	popB_ind7                	popB_ind5                	0.13062814343793533      
DIST(0)          	135          	16                 	15                 	popB_ind7                	popB_ind6                	0.14559357684436564      
DIST(0)          	136          	17                 	0                  	popB_ind8                	popA_ind1                	0.16567511205365193      
DIST(0)          	137          	17                 	1                  	popB_ind8                	popA_ind2                	0.12640773566882468      
DIST(0)          	138          	17                 	2                  	popB_ind8                	popA_ind3                	0.1708866189778657       
DIST(0)          	139          	17                 	3                  	popB_ind8                	popA_ind4                	0.23133671988418389      
DIST(0)          	140          	17                 	4                  	popB_ind8                	popA_ind5                	0.19225338752473922      
DIST(0)          	141          	17                 	5                  	popB_ind8                	popA_ind6                	0.19726913409919919      
DIST(0)          	142          	17                 	6                  	popB_ind8                	popA_ind7                	0.1479487393021007       
DIST(0)          	143          	17                 	7                  	popB_ind8                	popA_ind8                	0.17227982719132304      
DIST(0)          	144          	17                 	8                  	popB_ind8                	popA_ind9                	0.14959957714123032      
DIST(0)          	145          	17                 	9                  	popB_ind8                	popA_ind10               	0.15249983598143485      
DIST(0)          	146          	17                 	10                 	popB_ind8                	popB_ind1                	0.096185128202881956     
DIST(0)          	147          	17                 	11                 	popB_ind8                	popB_ind2                	0.12356963919485955      
DIST(0)          	148          	17                 	12                 	popB_ind8                	popB_ind3                	0.097396995678103115     
DIST(0)          	149          	17                 	13                 	popB_ind8                	popB_ind4                	0.1170702575416423       
DIST(0)          	150          	17                 	14                 	popB_ind8                	popB_ind5                	0.077202720953757881     
DIST(0)          	151          	17                 	15                 	popB_ind8                	popB_ind6                	0.12856793997731342      
DIST(0)          	152          	17                 	16                 	popB_ind8                	popB_ind7                	0.12011526516704371      
DIST(0)          	153          	18                 	0                  	popB_ind9                	popA_ind1                	0.12538366546380877      
DIST(0)          	154          	18                 	1                  	popB_ind9                	popA_ind2                	0.1536787031523098       
DIST(0)          	155          	18                 	2                  	popB_ind9                	popA_ind3                	0.10515210094250803      
DIST(0)          	156          	18                 	3                  	popB_ind9                	popA_ind4                	0.15781540012958445      
DIST(0)          	157          	18                 	4                  	popB_ind9                	popA_ind5                	0.12211001480676476      
DIST(0)          	158          	18                 	5                  	popB_ind9                	popA_ind6                	0.13742681971106568      
DIST(0)          	159          	18                 	6                  	popB_ind9                	popA_ind7                	0.11147969497017542      
DIST(0)          	160          	18                 	7                  	popB_ind9                	popA_ind8                	0.14821969503273411      
DIST(0)          	161          	18                 	8                  	popB_ind9                	popA_ind9                	0.095865949319817401     
DIST(0)          	162          	18                 	9                  	popB_ind9                	popA_ind10               	0.15796362529584287      
DIST(0)          	163          	18                 	10                 	popB_ind9                	popB_ind1                	0.11159237196376758      
DIST(0)          	164          	18                 	11                 	popB_ind9                	popB_ind2                	0.10566166338278885      
DIST(0)          	165          	18                 	12                 	popB_ind9                	popB_ind3                	0.1183479934870411       
DIST(0)          	166          	18                 	13                 	popB_ind9                	popB_ind4                	0.097778578496353533     
DIST(0)          	167          	18                 	14                 	popB_ind9                	popB_ind5                	0.092115999449317809     
DIST(0)          	168          	18                 	15                 	popB_ind9                	popB_ind6                	0.11123602906816318      
DIST(0)          	169          	18                 	16                 	popB_ind9                	popB_ind7                	0.079927548958855121     
DIST(0)          	170          	18                 	17                 	popB_ind9                	popB_ind8                	0.097300091716630593     
DIST(0)          	171          	19                 	0                  	popB_ind10               	popA_ind1                	0.2130766619311581       
DIST(0)          	172          	19                 	1                  	popB_ind10               	popA_ind2                	0.17037706765265032      
DIST(0)          	173          	19                 	2                  	popB_ind10               	popA_ind3                	0.14639700483239165      
DIST(0)          	174          	19                 	3                  	popB_ind10               	popA_ind4                	0.19214745347349027      
DIST(0)          	175          	19                 	4                  	popB_ind10               	popA_ind5                	0.13061699069084559      
DIST(0)          	176          	19                 	5                  	popB_ind10               	popA_ind6                	0.15657229526571761      
DIST(0)          	177          	19                 	6                  	popB_ind10               	popA_ind7                	0.1698556281188352       
DIST(0)          	178          	19                 	7                  	popB_ind10               	popA_ind8                	0.18171893336639938      
DIST(0)          	179          	19                 	8                  	popB_ind10               	popA_ind9                	0.16719664771957576      
DIST(0)          	180          	19                 	9                  	popB_ind10               	popA_ind10               	0.1634737129064138       
DIST(0)          	181          	19                 	10                 	popB_ind10               	popB_ind1                	0.12687343944407434      
DIST(0)          	182          	19                 	11                 	popB_ind10               	popB_ind2                	0.1443034076029559       
DIST(0)          	183          	19                 	12                 	popB_ind10               	popB_ind3                	0.071379753531828224     
DIST(0)          	184          	19                 	13                 	popB_ind10               	popB_ind4                	0.074913241182063398     
DIST(0)          	185          	19                 	14                 	popB_ind10               	popB_ind5                	0.11966181426920425      
DIST(0)          	186          	19                 	15                 	popB_ind10               	popB_ind6                	0.14635490020317696      
DIST(0)          	187          	19                 	16                 	popB_ind10               	popB_ind7                	0.093136998894170486     
DIST(0)          	188          	19                 	17                 	popB_ind10               	popB_ind8                	0.079366584664830098     
DIST(0)          	189          	19                 	18                 	popB_ind10               	popB_ind9                	0.1142753889662504       
DIST(0)          	190          	20                 	0                  	popC_ind1                	popA_ind1                	0.2015332824908172       
DIST(0)          	191          	20                 	1                  	popC_ind1                	popA_ind2                	0.21319717169203276      
DIST(0)          	192          	20                 	2                  	popC_ind1                	popA_ind3                	0.1993509842944885       
DIST(0)          	193          	20                 	3                  	popC_ind1                	popA_ind4                	0.22274743747405046      
DIST(0)          	194          	20                 	4                  	popC_ind1                	popA_ind5                	0.13482631394144037      
DIST(0)          	195          	20                 	5                  	popC_ind1                	popA_ind6                	0.22349667533920586      
DIST(0)          	196          	20                 	6                  	popC_ind1                	popA_ind7                	0.21226372044518818      
DIST(0)          	197          	20                 	7                  	popC_ind1                	popA_ind8                	0.15416149791707109      
DIST(0)          	198          	20                 	8                  	popC_ind1                	popA_ind9                	0.14624866514036186      
DIST(0)          	199          	20                 	9                  	popC_ind1                	popA_ind10               	0.19273603119366231      
DIST(0)          	200          	20                 	10                 	popC_ind1                	popB_ind1                	0.18992584418755429      
DIST(0)          	201          	20                 	11                 	popC_ind1                	popB_ind2                	0.13916576111370943      
DIST(0)          	202          	20                 	12                 	popC_ind1                	popB_ind3                	0.218718029433946        
DIST(0)          	203          	20                 	13                 	popC_ind1                	popB_ind4                	0.23894372190687277      
DIST(0)          	204          	20                 	14                 	popC_ind1                	popB_ind5                	0.17615682862512977      
DIST(0)          	205          	20                 	15                 	popC_ind1                	popB_ind6                	0.21229738138034671      
DIST(0)          	206          	20                 	16                 	popC_ind1                	popB_ind7                	0.25276437247717165      
DIST(0)          	207          	20                 	17                 	popC_ind1                	popB_ind8                	0.22832392155407685      
DIST(0)          	208          	20                 	18                 	popC_ind1                	popB_ind9                	0.20775044417262856      
DIST(0)          	209          	20                 	19                 	popC_ind1                	popB_ind10               	0.1441944134983002       
DIST(0)          	210          	21                 	0                  	popC_ind2                	popA_ind1                	0.22010928878735495      
DIST(0)          	211          	21                 	1                  	popC_ind2                	popA_ind2                	0.19596048314884856      
DIST(0)          	212          	21                 	2                  	popC_ind2                	popA_ind3                	0.13348044864677658      
DIST(0)          	213          	21                 	3                  	popC_ind2                	popA_ind4                	0.22768146183832755      
DIST(0)          	214          	21                 	4                  	popC_ind2                	popA_ind5                	0.17210693161572244      
DIST(0)          	215          	21                 	5                  	popC_ind2                	popA_ind6                	0.22974283767662557      
DIST(0)          	216          	21                 	6                  	popC_ind2                	popA_ind7                	0.20868695076057758      
DIST(0)          	217          	21                 	7                  	popC_ind2                	popA_ind8                	0.21980391670649729      
DIST(0)          	218          	21                 	8                  	popC_ind2                	popA_ind9                	0.20436847500812053      
DIST(0)          	219          	21                 	9                  	popC_ind2                	popA_ind10               	0.23773925327304762      
DIST(0)          	220          	21                 	10                 	popC_ind2                	popB_ind1                	0.2411019389417306       
DIST(0)          	221          	21                 	11                 	popC_ind2                	popB_ind2                	0.18274924278785148      
DIST(0)          	222          	21                 	12                 	popC_ind2                	popB_ind3                	0.22482126011468145      
DIST(0)          	223          	21                 	13                 	popC_ind2                	popB_ind4                	0.21303988719470016      
DIST(0)          	224          	21                 	14                 	popC_ind2                	popB_ind5                	0.24704618511131388      
DIST(0)          	225          	21                 	15                 	popC_ind2                	popB_ind6                	0.17241113608750086      
DIST(0)          	226          	21                 	16                 	popC_ind2                	popB_ind7                	0.22743761668678628      
DIST(0)          	227          	21                 	17                 	popC_ind2                	popB_ind8                	0.24500901223963262      
DIST(0)          	228          	21                 	18                 	popC_ind2                	popB_ind9                	0.17401907242445552      
DIST(0)          	229          	21                 	19                 	popC_ind2                	popB_ind10               	0.17641728721519201      
DIST(0)          	230          	21                 	20                 	popC_ind2                	popC_ind1                	0.085929952914894553     
DIST(0)          	231          	22                 	0                  	popC_ind3                	popA_ind1                	0.20992072881488485      
DIST(0)          	232          	22                 	1                  	popC_ind3                	popA_ind2                	0.19856025357945112      
DIST(0)          	233          	22                 	2                  	popC_ind3                	popA_ind3                	0.18137186054747748      
DIST(0)          	234          	22                 	3                  	popC_ind3                	popA_ind4                	0.22563648648046836      
DIST(0)          	235          	22                 	4                  	popC_ind3                	popA_ind5                	0.19355061640016907      
DIST(0)          	236          	22                 	5                  	popC_ind3                	popA_ind6                	0.20568966173104586      
DIST(0)          	237          	22                 	6                  	popC_ind3                	popA_ind7                	0.18097502863792983      
DIST(0)          	238          	22                 	7                  	popC_ind3                	popA_ind8                	0.19008316661633079      
DIST(0)          	239          	22                 	8                  	popC_ind3                	popA_ind9                	0.20914506792198181      
DIST(0)          	240          	22                 	9                  	popC_ind3                	popA_ind10               	0.196616965692434        
DIST(0)          	241          	22                 	10                 	popC_ind3                	popB_ind1                	0.20516226119894759      
DIST(0)          	242          	22                 	11                 	popC_ind3                	popB_ind2                	0.17419158370645471      
DIST(0)          	243          	22                 	12                 	popC_ind3                	popB_ind3                	0.22366785596444297      
DIST(0)          	244          	22                 	13                 	popC_ind3                	popB_ind4                	0.23600090782376576      
DIST(0)          	245          	22                 	14                 	popC_ind3                	popB_ind5                	0.24198609467351978      
DIST(0)          	246          	22                 	15                 	popC_ind3                	popB_ind6                	0.19888742525337738      
DIST(0)          	247          	22                 	16                 	popC_ind3                	popB_ind7                	0.21567635282876946      
DIST(0)          	248          	22                 	17                 	popC_ind3                	popB_ind8                	0.20019062125200141      
DIST(0)          	249          	22                 	18                 	popC_ind3                	popB_ind9                	0.20930096488445013      
DIST(0)          	250          	22                 	19                 	popC_ind3                	popB_ind10               	0.19961998820517368      
DIST(0)          	251          	22                 	20                 	popC_ind3                	popC_ind1                	0.064587240380522284     
DIST(0)          	252          	22                 	21                 	popC_ind3                	popC_ind2                	0.083556136378852702     
DIST(0)          	253          	23                 	0                  	popC_ind4                	popA_ind1                	0.18011973954276256      
DIST(0)          	254          	23                 	1                  	popC_ind4                	popA_ind2                	0.17799342432036944      
DIST(0)          	255          	23                 	2                  	popC_ind4                	popA_ind3                	0.14461749966416013      
DIST(0)          	256          	23                 	3                  	popC_ind4                	popA_ind4                	0.19707754359541055      
DIST(0)          	257          	23                 	4                  	popC_ind4                	popA_ind5                	0.081157311261563433     
DIST(0)          	258          	23                 	5                  	popC_ind4                	popA_ind6                	0.19920225108790207      
DIST(0)          	259          	23                 	6                  	popC_ind4                	popA_ind7                	0.23824256705359415      
DIST(0)          	260          	23                 	7                  	popC_ind4                	popA_ind8                	0.12737489517049935      
DIST(0)          	261          	23                 	8                  	popC_ind4                	popA_ind9                	0.13183984482082817      
DIST(0)          	262          	23                 	9                  	popC_ind4                	popA_ind10               	0.16902390438926701      
DIST(0)          	263          	23                 	10                 	popC_ind4                	popB_ind1                	0.16848696867103974      
DIST(0)          	264          	23                 	11                 	popC_ind4                	popB_ind2                	0.14855713001276011      
DIST(0)          	265          	23                 	12                 	popC_ind4                	popB_ind3                	0.22618623471090815      
DIST(0)          	266          	23                 	13                 	popC_ind4                	popB_ind4                	0.19733227421567026      
DIST(0)          	267          	23                 	14                 	popC_ind4                	popB_ind5                	0.219284944468256        
DIST(0)          	268          	23                 	15                 	popC_ind4                	popB_ind6                	0.14381946282788433      
DIST(0)          	269          	23                 	16                 	popC_ind4                	popB_ind7                	0.18400919938333199      
DIST(0)          	270          	23                 	17                 	popC_ind4                	popB_ind8                	0.21535393907604786      
DIST(0)          	271          	23                 	18                 	popC_ind4                	popB_ind9                	0.17428643282842965      
DIST(0)          	272          	23                 	19                 	popC_ind4                	popB_ind10               	0.16779312578307787      
DIST(0)          	273          	23                 	20                 	popC_ind4                	popC_ind1                	0.088996394857584024     
DIST(0)          	274          	23                 	21                 	popC_ind4                	popC_ind2                	0.063669140840762459     
DIST(0)          	275          	23                 	22                 	popC_ind4                	popC_ind3                	0.087365595598285273     
DIST(0)          	276          	24                 	0                  	popC_ind5                	popA_ind1                	0.20325323130980855      
DIST(0)          	277          	24                 	1                  	popC_ind5                	popA_ind2                	0.14702450642836398      
DIST(0)          	278          	24                 	2                  	popC_ind5                	popA_ind3                	0.12378919856248134      
DIST(0)          	279          	24                 	3                  	popC_ind5                	popA_ind4                	0.16490775905017863      
DIST(0)          	280          	24                 	4                  	popC_ind5                	popA_ind5                	0.12756045244004605      
DIST(0)          	281          	24                 	5                  	popC_ind5                	popA_ind6                	0.1527573957747721       
DIST(0)          	282          	24                 	6                  	popC_ind5                	popA_ind7                	0.14767875278433937      
DIST(0)          	283          	24                 	7                  	popC_ind5                	popA_ind8                	0.11134245564017317      
DIST(0)          	284          	24                 	8                  	popC_ind5                	popA_ind9                	0.15457974651569453      
DIST(0)          	285          	24                 	9                  	popC_ind5                	popA_ind10               	0.16103514198867447      
DIST(0)          	286          	24                 	10                 	popC_ind5                	popB_ind1                	0.18153955539439673      
DIST(0)          	287          	24                 	11                 	popC_ind5                	popB_ind2                	0.12645621735530055      
DIST(0)          	288          	24                 	12                 	popC_ind5                	popB_ind3                	0.21219604150214655      
DIST(0)          	289          	24                 	13                 	popC_ind5                	popB_ind4                	0.22434630118986329      
DIST(0)          	290          	24                 	14                 	popC_ind5                	popB_ind5                	0.21254986108129215      
DIST(0)          	291          	24                 	15                 	popC_ind5                	popB_ind6                	0.16881940483549152      
DIST(0)          	292          	24                 	16                 	popC_ind5                	popB_ind7                	0.14723056087689204      
DIST(0)          	293          	24                 	17                 	popC_ind5                	popB_ind8                	0.20735358341427793      
DIST(0)          	294          	24                 	18                 	popC_ind5                	popB_ind9                	0.15812908243462714      
DIST(0)          	295          	24                 	19                 	popC_ind5                	popB_ind10               	0.12767489028796095      
DIST(0)          	296          	24                 	20                 	popC_ind5                	popC_ind1                	0.038373039060623441     
DIST(0)          	297          	24                 	21                 	popC_ind5                	popC_ind2                	0.092059644589807504     
DIST(0)          	298          	24                 	22                 	popC_ind5                	popC_ind3                	0.074413818107818733     
DIST(0)          	299          	24                 	23                 	popC_ind5                	popC_ind4                	0.053726777929373533     
DIST(0)          	300          	25                 	0                  	popC_ind6                	popA_ind1                	0.20984031673621631      
DIST(0)          	301          	25                 	1                  	popC_ind6                	popA_ind2                	0.17442468605602207      
DIST(0)          	302          	25                 	2                  	popC_ind6                	popA_ind3                	0.2589571139309198       
DIST(0)          	303          	25                 	3                  	popC_ind6                	popA_ind4                	0.2116313684226776       
DIST(0)          	304          	25                 	4                  	popC_ind6                	popA_ind5                	0.16355731741209178      
DIST(0)          	305          	25                 	5                  	popC_ind6                	popA_ind6                	0.26326882256131295      
DIST(0)          	306          	25                 	6                  	popC_ind6                	popA_ind7                	0.26393090198139946      
DIST(0)          	307          	25                 	7                  	popC_ind6                	popA_ind8                	0.16815725826357403      
DIST(0)          	308          	25                 	8                  	popC_ind6                	popA_ind9                	0.22714663644873606      
DIST(0)          	309          	25                 	9                  	popC_ind6                	popA_ind10               	0.24610165008253948      
DIST(0)          	310          	25                 	10                 	popC_ind6                	popB_ind1                	0.24285448097296569      
DIST(0)          	311          	25                 	11                 	popC_ind6                	popB_ind2                	0.18794323710484434      
DIST(0)          	312          	25                 	12                 	popC_ind6                	popB_ind3                	0.24428109312872376      
DIST(0)          	313          	25                 	13                 	popC_ind6                	popB_ind4                	0.18087135708310498      
DIST(0)          	314          	25                 	14                 	popC_ind6                	popB_ind5                	0.18653439343397094      
DIST(0)          	315          	25                 	15                 	popC_ind6                	popB_ind6                	0.20582109677591914      
DIST(0)          	316          	25                 	16                 	popC_ind6                	popB_ind7                	0.22320203943759706      
DIST(0)          	317          	25                 	17                 	popC_ind6                	popB_ind8                	0.21083984955048846      
DIST(0)          	318          	25                 	18                 	popC_ind6                	popB_ind9                	0.18504443247286206      
DIST(0)          	319          	25                 	19                 	popC_ind6                	popB_ind10               	0.1909788577390174       
DIST(0)          	320          	25                 	20                 	popC_ind6                	popC_ind1                	0.061336527554696936     
DIST(0)          	321          	25                 	21                 	popC_ind6                	popC_ind2                	0.060911166662707003     
DIST(0)          	322          	25                 	22                 	popC_ind6                	popC_ind3                	0.088665456548880486     
DIST(0)          	323          	25                 	23                 	popC_ind6                	popC_ind4                	0.074081134186570141     
DIST(0)          	324          	25                 	24                 	popC_ind6                	popC_ind5                	0.053431671951617909     
DIST(0)          	325          	26                 	0                  	popC_ind7                	popA_ind1                	0.23655857744765843      
DIST(0)          	326          	26                 	1                  	popC_ind7                	popA_ind2                	0.24830740085268679      
DIST(0)          	327          	26                 	2                  	popC_ind7                	popA_ind3                	0.2216585993132843       
DIST(0)          	328          	26                 	3                  	popC_ind7                	popA_ind4                	0.18967820697965593      
DIST(0)          	329          	26                 	4                  	popC_ind7                	popA_ind5                	0.2013899855208609       
DIST(0)          	330          	26                 	5                  	popC_ind7                	popA_ind6                	0.28795560981523416      
DIST(0)          	331          	26                 	6                  	popC_ind7                	popA_ind7                	0.26873541148604974      
DIST(0)          	332          	26                 	7                  	popC_ind7                	popA_ind8                	0.19675588091388391      
DIST(0)          	333          	26                 	8                  	popC_ind7                	popA_ind9                	0.22811461996768145      
DIST(0)          	334          	26                 	9                  	popC_ind7                	popA_ind10               	0.20867271024605338      
DIST(0)          	335          	26                 	10                 	popC_ind7                	popB_ind1                	0.24061109001033371      
DIST(0)          	336          	26                 	11                 	popC_ind7                	popB_ind2                	0.19877923217784799      
DIST(0)          	337          	26                 	12                 	popC_ind7                	popB_ind3                	0.21547008840665582      
DIST(0)          	338          	26                 	13                 	popC_ind7                	popB_ind4                	0.19767769172169311      
DIST(0)          	339          	26                 	14                 	popC_ind7                	popB_ind5                	0.27509638561703642      
DIST(0)          	340          	26                 	15                 	popC_ind7                	popB_ind6                	0.25374463338126252      
DIST(0)          	341          	26                 	16                 	popC_ind7                	popB_ind7                	0.17857892695742406      
DIST(0)          	342          	26                 	17                 	popC_ind7                	popB_ind8                	0.22978234231581426      
DIST(0)          	343          	26                 	18                 	popC_ind7                	popB_ind9                	0.22425103724399484      
DIST(0)          	344          	26                 	19                 	popC_ind7                	popB_ind10               	0.19836430073578473      
DIST(0)          	345          	26                 	20                 	popC_ind7                	popC_ind1                	0.080468945192271019     
DIST(0)          	346          	26                 	21                 	popC_ind7                	popC_ind2                	0.077573278036358176     
DIST(0)          	347          	26                 	22                 	popC_ind7                	popC_ind3                	0.096624243274473623     
DIST(0)          	348          	26                 	23                 	popC_ind7                	popC_ind4                	0.078420141511245114     
DIST(0)          	349          	26                 	24                 	popC_ind7                	popC_ind5                	0.085692851086862684     
DIST(0)          	350          	26                 	25                 	popC_ind7                	popC_ind6                	0.094015717217159234     
DIST(0)          	351          	27                 	0                  	popC_ind8                	popA_ind1                	0.23171725978169411      
DIST(0)          	352          	27                 	1                  	popC_ind8                	popA_ind2                	0.22090372507239436      
DIST(0)          	353          	27                 	2                  	popC_ind8                	popA_ind3                	0.20809777426540216      
DIST(0)          	354          	27                 	3                  	popC_ind8                	popA_ind4                	0.19897184497235865      
DIST(0)          	355          	27                 	4                  	popC_ind8                	popA_ind5                	0.2098146614756794       
DIST(0)          	356          	27                 	5                  	popC_ind8                	popA_ind6                	0.23941607366271042      
DIST(0)          	357          	27                 	6                  	popC_ind8                	popA_ind7                	0.20018477116636771      
DIST(0)          	358          	27                 	7                  	popC_ind8                	popA_ind8                	0.17347278125565394      
DIST(0)          	359          	27                 	8                  	popC_ind8                	popA_ind9                	0.20419846762542376      
DIST(0)          	360          	27                 	9                  	popC_ind8                	popA_ind10               	0.21781060486206261      
DIST(0)          	361          	27                 	10                 	popC_ind8                	popB_ind1                	0.22480874276002177      
DIST(0)          	362          	27                 	11                 	popC_ind8                	popB_ind2                	0.27709420638695614      
DIST(0)          	363          	27                 	12                 	popC_ind8                	popB_ind3                	0.22787476566118275      
DIST(0)          	364          	27                 	13                 	popC_ind8                	popB_ind4                	0.19083714046180256      
DIST(0)          	365          	27                 	14                 	popC_ind8                	popB_ind5                	0.30645922919112667      
DIST(0)          	366          	27                 	15                 	popC_ind8                	popB_ind6                	0.23059273142635181      
DIST(0)          	367          	27                 	16                 	popC_ind8                	popB_ind7                	0.23069574344760835      
DIST(0)          	368          	27                 	17                 	popC_ind8                	popB_ind8                	0.26584523918957537      
DIST(0)          	369          	27                 	18                 	popC_ind8                	popB_ind9                	0.22086612680901901      
DIST(0)          	370          	27                 	19                 	popC_ind8                	popB_ind10               	0.19970167979323389      
DIST(0)          	371          	27                 	20                 	popC_ind8                	popC_ind1                	0.038367044363519641     
DIST(0)          	372          	27                 	21                 	popC_ind8                	popC_ind2                	0.057603914036350262     
DIST(0)          	373          	27                 	22                 	popC_ind8                	popC_ind3                	0.10002365391216407      
DIST(0)          	374          	27                 	23                 	popC_ind8                	popC_ind4                	0.074811407190030099     
DIST(0)          	375          	27                 	24                 	popC_ind8                	popC_ind5                	0.049541727827786639     
DIST(0)          	376          	27                 	25                 	popC_ind8                	popC_ind6                	0.083376674419933378     
DIST(0)          	377          	27                 	26                 	popC_ind8                	popC_ind7                	0.10231720437086898      
DIST(0)          	378          	28                 	0                  	popC_ind9                	popA_ind1                	0.27893889901517865      
DIST(0)          	379          	28                 	1                  	popC_ind9                	popA_ind2                	0.21410309795138521      
DIST(0)          	380          	28                 	2                  	popC_ind9                	popA_ind3                	0.19286296697013489      
DIST(0)          	381          	28                 	3                  	popC_ind9                	popA_ind4                	0.18686749308825451      
DIST(0)          	382          	28                 	4                  	popC_ind9                	popA_ind5                	0.15456903913146047      
DIST(0)          	383          	28                 	5                  	popC_ind9                	popA_ind6                	0.21061187153154842      
DIST(0)          	384          	28                 	6                  	popC_ind9                	popA_ind7                	0.22837868173009898      
DIST(0)          	385          	28                 	7                  	popC_ind9                	popA_ind8                	0.1637846029129669       
DIST(0)          	386          	28                 	8                  	popC_ind9                	popA_ind9                	0.17143383683045679      
DIST(0)          	387          	28                 	9                  	popC_ind9                	popA_ind10               	0.14333689465982716      
DIST(0)          	388          	28                 	10                 	popC_ind9                	popB_ind1                	0.2403157110690059       
DIST(0)          	389          	28                 	11                 	popC_ind9                	popB_ind2                	0.19753342149877423      
DIST(0)          	390          	28                 	12                 	popC_ind9                	popB_ind3                	0.24909517760238462      
DIST(0)          	391          	28                 	13                 	popC_ind9                	popB_ind4                	0.21178624032341026      
DIST(0)          	392          	28                 	14                 	popC_ind9                	popB_ind5                	0.23701875375463244      
DIST(0)          	393          	28                 	15                 	popC_ind9                	popB_ind6                	0.20842888691417302      
DIST(0)          	394          	28                 	16                 	popC_ind9                	popB_ind7                	0.27955586354207185      
DIST(0)          	395          	28                 	17                 	popC_ind9                	popB_ind8                	0.29242445900054809      
DIST(0)          	396          	28                 	18                 	popC_ind9                	popB_ind9                	0.25257668494528768      
DIST(0)          	397          	28                 	19                 	popC_ind9                	popB_ind10               	0.19507042811246808      
DIST(0)          	398          	28                 	20                 	popC_ind9                	popC_ind1                	0.062541815271601547     
DIST(0)          	399          	28                 	21                 	popC_ind9                	popC_ind2                	0.096198048711564468     
DIST(0)          	400          	28                 	22                 	popC_ind9                	popC_ind3                	0.1280191774205563       
DIST(0)          	401          	28                 	23                 	popC_ind9                	popC_ind4                	0.095481720882076723     
DIST(0)          	402          	28                 	24                 	popC_ind9                	popC_ind5                	0.046700724870551005     
DIST(0)          	403          	28                 	25                 	popC_ind9                	popC_ind6                	0.078023295539278489     
DIST(0)          	404          	28                 	26                 	popC_ind9                	popC_ind7                	0.094692791475228327     
DIST(0)          	405          	28                 	27                 	popC_ind9                	popC_ind8                	0.094952695531240541     
DIST(0)          	406          	29                 	0                  	popC_ind10               	popA_ind1                	0.26359692965102521      
DIST(0)          	407          	29                 	1                  	popC_ind10               	popA_ind2                	0.24990291569905126      
DIST(0)          	408          	29                 	2                  	popC_ind10               	popA_ind3                	0.24506414677891708      
DIST(0)          	409          	29                 	3                  	popC_ind10               	popA_ind4                	0.27575413391200859      
DIST(0)          	410          	29                 	4                  	popC_ind10               	popA_ind5                	0.24426070118934098      
DIST(0)          	411          	29                 	5                  	popC_ind10               	popA_ind6                	0.22138250336297252      
DIST(0)          	412          	29                 	6                  	popC_ind10               	popA_ind7                	0.22809973495139552      
DIST(0)          	413          	29                 	7                  	popC_ind10               	popA_ind8                	0.22615354130733101      
DIST(0)          	414          	29                 	8                  	popC_ind10               	popA_ind9                	0.24653119815044355      
DIST(0)          	415          	29                 	9                  	popC_ind10               	popA_ind10               	0.2294934646127805       
DIST(0)          	416          	29                 	10                 	popC_ind10               	popB_ind1                	0.20338124312427919      
DIST(0)          	417          	29                 	11                 	popC_ind10               	popB_ind2                	0.23936921589388946      
DIST(0)          	418          	29                 	12                 	popC_ind10               	popB_ind3                	0.24991241331889907      
DIST(0)          	419          	29                 	13                 	popC_ind10               	popB_ind4                	0.22499187380586272      
DIST(0)          	420          	29                 	14                 	popC_ind10               	popB_ind5                	0.27313820987458215      
DIST(0)          	421          	29                 	15                 	popC_ind10               	popB_ind6                	0.25544194022353872      
DIST(0)          	422          	29                 	16                 	popC_ind10               	popB_ind7                	0.2754649220079905       
DIST(0)          	423          	29                 	17                 	popC_ind10               	popB_ind8                	0.2684169366788503       
DIST(0)          	424          	29                 	18                 	popC_ind10               	popB_ind9                	0.21351601714660881      
DIST(0)          	425          	29                 	19                 	popC_ind10               	popB_ind10               	0.19710585973442807      
DIST(0)          	426          	29                 	20                 	popC_ind10               	popC_ind1                	0.084133763987661189     
DIST(0)          	427          	29                 	21                 	popC_ind10               	popC_ind2                	0.089002603093127908     
DIST(0)          	428          	29                 	22                 	popC_ind10               	popC_ind3                	0.090197021640572694     
DIST(0)          	429          	29                 	23                 	popC_ind10               	popC_ind4                	0.070476984843124091     
DIST(0)          	430          	29                 	24                 	popC_ind10               	popC_ind5                	0.069405838498509606     
DIST(0)          	431          	29                 	25                 	popC_ind10               	popC_ind6                	0.081919128661087498     
DIST(0)          	432          	29                 	26                 	popC_ind10               	popC_ind7                	0.10398249075452988      
DIST(0)          	433          	29                 	27                 	popC_ind10               	popC_ind8                	0.088913855721502233     
DIST(0)          	434          	29                 	28                 	popC_ind10               	popC_ind9                	0.061781873583093873     
DIST(0)          	435          	30                 	0                  	popD_ind1                	popA_ind1                	0.1557687437607343       
DIST(0)          	436          	30                 	1                  	popD_ind1                	popA_ind2                	0.18354603047325468      
DIST(0)          	437          	30                 	2                  	popD_ind1                	popA_ind3                	0.15548244825433166      
DIST(0)          	438          	30                 	3                  	popD_ind1                	popA_ind4                	0.23636189166322444      
DIST(0)          	439          	30                 	4                  	popD_ind1                	popA_ind5                	0.16065407838892959      
DIST(0)          	440          	30                 	5                  	popD_ind1                	popA_ind6                	0.25146596833430723      
DIST(0)          	441          	30                 	6                  	popD_ind1                	popA_ind7                	0.18426987979109452      
DIST(0)          	442          	30                 	7                  	popD_ind1                	popA_ind8                	0.1424256545264975       
DIST(0)          	443          	30                 	8                  	popD_ind1                	popA_ind9                	0.18239425024363842      
DIST(0)          	444          	30                 	9                  	popD_ind1                	popA_ind10               	0.20453525675933962      
DIST(0)          	445          	30                 	10                 	popD_ind1                	popB_ind1                	0.20477225865304255      
DIST(0)          	446          	30                 	11                 	popD_ind1                	popB_ind2                	0.20888693041062814      
DIST(0)          	447          	30                 	12                 	popD_ind1                	popB_ind3                	0.19351193017586943      
DIST(0)          	448          	30                 	13                 	popD_ind1                	popB_ind4                	0.199913543115711        
DIST(0)          	449          	30                 	14                 	popD_ind1                	popB_ind5                	0.24126918206266978      
DIST(0)          	450          	30                 	15                 	popD_ind1                	popB_ind6                	0.17394113001215428      
DIST(0)          	451          	30                 	16                 	popD_ind1                	popB_ind7                	0.18805765501821292      
DIST(0)          	452          	30                 	17                 	popD_ind1                	popB_ind8                	0.18541141583156168      
DIST(0)          	453          	30                 	18                 	popD_ind1                	popB_ind9                	0.15971641445464083      
DIST(0)          	454          	30                 	19                 	popD_ind1                	popB_ind10               	0.20158372573305072      
DIST(0)          	455          	30                 	20                 	popD_ind1                	popC_ind1                	0.12984999798110172      
DIST(0)          	456          	30                 	21                 	popD_ind1                	popC_ind2                	0.093924844356783055     
DIST(0)          	457          	30                 	22                 	popD_ind1                	popC_ind3                	0.14958527171031463      
DIST(0)          	458          	30                 	23                 	popD_ind1                	popC_ind4                	0.087817162380279215     
DIST(0)          	459          	30                 	24                 	popD_ind1                	popC_ind5                	0.074551563154910017     
DIST(0)          	460          	30                 	25                 	popD_ind1                	popC_ind6                	0.11993332403573877      
DIST(0)          	461          	30                 	26                 	popD_ind1                	popC_ind7                	0.16517261220996154      
DIST(0)          	462          	30                 	27                 	popD_ind1                	popC_ind8                	0.11194600643558583      
DIST(0)          	463          	30                 	28                 	popD_ind1                	popC_ind9                	0.13796409469097348      
DIST(0)          	464          	30                 	29                 	popD_ind1                	popC_ind10               	0.16068347267613489      
DIST(0)          	465          	31                 	0                  	popD_ind2                	popA_ind1                	0.2252919853265799       
DIST(0)          	466          	31                 	1                  	popD_ind2                	popA_ind2                	0.15273144787650705      
DIST(0)          	467          	31                 	2                  	popD_ind2                	popA_ind3                	0.23137037387074691      
DIST(0)          	468          	31                 	3                  	popD_ind2                	popA_ind4                	0.21336392171997742      
DIST(0)          	469          	31                 	4                  	popD_ind2                	popA_ind5                	0.18413506961360532      
DIST(0)          	470          	31                 	5                  	popD_ind2                	popA_ind6                	0.20093618363776045      
DIST(0)          	471          	31                 	6                  	popD_ind2                	popA_ind7                	0.17424317963019315      
DIST(0)          	472          	31                 	7                  	popD_ind2                	popA_ind8                	0.22066812525444179      
DIST(0)          	473          	31                 	8                  	popD_ind2                	popA_ind9                	0.16343172198789713      
DIST(0)          	474          	31                 	9                  	popD_ind2                	popA_ind10               	0.2683787525760486       
DIST(0)          	475          	31                 	10                 	popD_ind2                	popB_ind1                	0.2176982159169541       
DIST(0)          	476          	31                 	11                 	popD_ind2                	popB_ind2                	0.22592371461100735      
DIST(0)          	477          	31                 	12                 	popD_ind2                	popB_ind3                	0.27272862947566057      
DIST(0)          	478          	31                 	13                 	popD_ind2                	popB_ind4                	0.20245036636524577      
DIST(0)          	479          	31                 	14                 	popD_ind2                	popB_ind5                	0.25020944057727357      
DIST(0)          	480          	31                 	15                 	popD_ind2                	popB_ind6                	0.26265371451431879      
DIST(0)          	481          	31                 	16                 	popD_ind2                	popB_ind7                	0.2265333001434483       
DIST(0)          	482          	31                 	17                 	popD_ind2                	popB_ind8                	0.23948524980857971      
DIST(0)          	483          	31                 	18                 	popD_ind2                	popB_ind9                	0.21237658573738638      
DIST(0)          	484          	31                 	19                 	popD_ind2                	popB_ind10               	0.19899818404453437      
DIST(0)          	485          	31                 	20                 	popD_ind2                	popC_ind1                	0.16846483892402273      
DIST(0)          	486          	31                 	21                 	popD_ind2                	popC_ind2                	0.166728731623842        
DIST(0)          	487          	31                 	22                 	popD_ind2                	popC_ind3                	0.16175438219767974      
DIST(0)          	488          	31                 	23                 	popD_ind2                	popC_ind4                	0.082048556450437959     
DIST(0)          	489          	31                 	24                 	popD_ind2                	popC_ind5                	0.065738464466928354     
DIST(0)          	490          	31                 	25                 	popD_ind2                	popC_ind6                	0.12302478124949576      
DIST(0)          	491          	31                 	26                 	popD_ind2                	popC_ind7                	0.1624192388846194       
DIST(0)          	492          	31                 	27                 	popD_ind2                	popC_ind8                	0.17524499344050826      
DIST(0)          	493          	31                 	28                 	popD_ind2                	popC_ind9                	0.11883894416167089      
DIST(0)          	494          	31                 	29                 	popD_ind2                	popC_ind10               	0.17443339134304406      
DIST(0)          	495          	31                 	30                 	popD_ind2                	popD_ind1                	0.064576161974153279     
DIST(0)          	496          	32                 	0                  	popD_ind3                	popA_ind1                	0.22101561676569037      
DIST(0)          	497          	32                 	1                  	popD_ind3                	popA_ind2                	0.091315068731080598     
DIST(0)          	498          	32                 	2                  	popD_ind3                	popA_ind3                	0.17047831749440684      
DIST(0)          	499          	32                 	3                  	popD_ind3                	popA_ind4                	0.18281749583573956      
DIST(0)          	500          	32                 	4                  	popD_ind3                	popA_ind5                	0.19268107311000748      
DIST(0)          	501          	32                 	5                  	popD_ind3                	popA_ind6                	0.22939261178266604      
DIST(0)          	502          	32                 	6                  	popD_ind3                	popA_ind7                	0.17338542867880857      
DIST(0)          	503          	32                 	7                  	popD_ind3                	popA_ind8                	0.19837001276276728      
DIST(0)          	504          	32                 	8                  	popD_ind3                	popA_ind9                	0.15008270114662817      
DIST(0)          	505          	32                 	9                  	popD_ind3                	popA_ind10               	0.18563201568799226      
DIST(0)          	506          	32                 	10                 	popD_ind3                	popB_ind1                	0.21893510830319929      
DIST(0)          	507          	32                 	11                 	popD_ind3                	popB_ind2                	0.18349102271350159      
DIST(0)          	508          	32                 	12                 	popD_ind3                	popB_ind3                	0.21277993095707076      
DIST(0)          	509          	32                 	13                 	popD_ind3                	popB_ind4                	0.17493860817319762      
DIST(0)          	510          	32                 	14                 	popD_ind3                	popB_ind5                	0.19796806509855142      
DIST(0)          	511          	32                 	15                 	popD_ind3                	popB_ind6                	0.17308135954369908      
DIST(0)          	512          	32                 	16                 	popD_ind3                	popB_ind7                	0.20618120362975073      
DIST(0)          	513          	32                 	17                 	popD_ind3                	popB_ind8                	0.1463432015373973       
DIST(0)          	514          	32                 	18                 	popD_ind3                	popB_ind9                	0.17480534081730839      
DIST(0)          	515          	32                 	19                 	popD_ind3                	popB_ind10               	0.17374194623019329      
DIST(0)          	516          	32                 	20                 	popD_ind3                	popC_ind1                	0.10930288882255604      
DIST(0)          	517          	32                 	21                 	popD_ind3                	popC_ind2                	0.10846395089436804      
DIST(0)          	518          	32                 	22                 	popD_ind3                	popC_ind3                	0.12783401100466404      
DIST(0)          	519          	32                 	23                 	popD_ind3                	popC_ind4                	0.11469374835607198      
DIST(0)          	520          	32                 	24                 	popD_ind3                	popC_ind5                	0.081969681336522646     
DIST(0)          	521          	32                 	25                 	popD_ind3                	popC_ind6                	0.13636551983836287      
DIST(0)          	522          	32                 	26                 	popD_ind3                	popC_ind7                	0.16349176700277329      
DIST(0)          	523          	32                 	27                 	popD_ind3                	popC_ind8                	0.11927721961442059      
DIST(0)          	524          	32                 	28                 	popD_ind3                	popC_ind9                	0.12640741666249491      
DIST(0)          	525          	32                 	29                 	popD_ind3                	popC_ind10               	0.20027109942585425      
DIST(0)          	526          	32                 	30                 	popD_ind3                	popD_ind1                	0.055549174616355043     
DIST(0)          	527          	32                 	31                 	popD_ind3                	popD_ind2                	0.084741591710787309     
DIST(0)          	528          	33                 	0                  	popD_ind4                	popA_ind1                	0.20867651590808914      
DIST(0)          	529          	33                 	1                  	popD_ind4                	popA_ind2                	0.1586076732132409       
DIST(0)          	530          	33                 	2                  	popD_ind4                	popA_ind3                	0.20642891361451218      
DIST(0)          	531          	33                 	3                  	popD_ind4                	popA_ind4                	0.19888747055935219      
DIST(0)          	532          	33                 	4                  	popD_ind4                	popA_ind5                	0.18392410152495564      
DIST(0)          	533          	33                 	5                  	popD_ind4                	popA_ind6                	0.22144594242960736      
DIST(0)          	534          	33                 	6                  	popD_ind4                	popA_ind7                	0.17902726726923524      
DIST(0)          	535          	33                 	7                  	popD_ind4                	popA_ind8                	0.20990936913047709      
DIST(0)          	536          	33                 	8                  	popD_ind4                	popA_ind9                	0.19507804468020856      
DIST(0)          	537          	33                 	9                  	popD_ind4                	popA_ind10               	0.18038684264531946      
DIST(0)          	538          	33                 	10                 	popD_ind4                	popB_ind1                	0.23919289155923101      
DIST(0)          	539          	33                 	11                 	popD_ind4                	popB_ind2                	0.15478947561292911      
DIST(0)          	540          	33                 	12                 	popD_ind4                	popB_ind3                	0.19941233218382945      
DIST(0)          	541          	33                 	13                 	popD_ind4                	popB_ind4                	0.20573764340109796      
DIST(0)          	542          	33                 	14                 	popD_ind4                	popB_ind5                	0.25391447973445547      
DIST(0)          	543          	33                 	15                 	popD_ind4                	popB_ind6                	0.23530740789842947      
DIST(0)          	544          	33                 	16                 	popD_ind4                	popB_ind7                	0.22514748421643299      
DIST(0)          	545          	33                 	17                 	popD_ind4                	popB_ind8                	0.21471598894825331      
DIST(0)          	546          	33                 	18                 	popD_ind4                	popB_ind9                	0.19174464449363682      
DIST(0)          	547          	33                 	19                 	popD_ind4                	popB_ind10               	0.21651897983173204      
DIST(0)          	548          	33                 	20                 	popD_ind4                	popC_ind1                	0.14684020519339983      
DIST(0)          	549          	33                 	21                 	popD_ind4                	popC_ind2                	0.12154685529474632      
DIST(0)          	550          	33                 	22                 	popD_ind4                	popC_ind3                	0.16840734321912113      
DIST(0)          	551          	33                 	23                 	popD_ind4                	popC_ind4                	0.1310283516519139       
DIST(0)          	552          	33                 	24                 	popD_ind4                	popC_ind5                	0.10528296756241991      
DIST(0)          	553          	33                 	25                 	popD_ind4                	popC_ind6                	0.092455188594602639     
DIST(0)          	554          	33                 	26                 	popD_ind4                	popC_ind7                	0.15531964686520139      
DIST(0)          	555          	33                 	27                 	popD_ind4                	popC_ind8                	0.17810137419091548      
DIST(0)          	556          	33                 	28                 	popD_ind4                	popC_ind9                	0.11642410183818176      
DIST(0)          	557          	33                 	29                 	popD_ind4                	popC_ind10               	0.16329440387638294      
DIST(0)          	558          	33                 	30                 	popD_ind4                	popD_ind1                	0.040291159308572445     
DIST(0)          	559          	33                 	31                 	popD_ind4                	popD_ind2                	0.071182101174123089     
DIST(0)          	560          	33                 	32                 	popD_ind4                	popD_ind3                	0.050850306031042618     
DIST(0)          	561          	34                 	0                  	popD_ind5                	popA_ind1                	0.2288501073715965       
DIST(0)          	562          	34                 	1                  	popD_ind5                	popA_ind2                	0.23106126738994098      
DIST(0)          	563          	34                 	2                  	popD_ind5                	popA_ind3                	0.21736566698621612      
DIST(0)          	564          	34                 	3                  	popD_ind5                	popA_ind4                	0.21130195979700117      
DIST(0)          	565          	34                 	4                  	popD_ind5                	popA_ind5                	0.25371929719601982      
DIST(0)          	566          	34                 	5                  	popD_ind5                	popA_ind6                	0.25158744905931074      
DIST(0)          	567          	34                 	6                  	popD_ind5                	popA_ind7                	0.23206171111174118      
DIST(0)          	568          	34                 	7                  	popD_ind5                	popA_ind8                	0.22433396965282579      
DIST(0)          	569          	34                 	8                  	popD_ind5                	popA_ind9                	0.20455029119506379      
DIST(0)          	570          	34                 	9                  	popD_ind5                	popA_ind10               	0.22525837630830783      
DIST(0)          	571          	34                 	10                 	popD_ind5                	popB_ind1                	0.22001406799952861      
DIST(0)          	572          	34                 	11                 	popD_ind5                	popB_ind2                	0.17962672021236251      
DIST(0)          	573          	34                 	12                 	popD_ind5                	popB_ind3                	0.19070710283874676      
DIST(0)          	574          	34                 	13                 	popD_ind5                	popB_ind4                	0.17960029739689223      
DIST(0)          	575          	34                 	14                 	popD_ind5                	popB_ind5                	0.24667554847277942      
DIST(0)          	576          	34                 	15                 	popD_ind5                	popB_ind6                	0.24720539059618932      
DIST(0)          	577          	34                 	16                 	popD_ind5                	popB_ind7                	0.22253210946423452      
DIST(0)          	578          	34                 	17                 	popD_ind5                	popB_ind8                	0.22024858207637088      
DIST(0)          	579          	34                 	18                 	popD_ind5                	popB_ind9                	0.22087491212067947      
DIST(0)          	580          	34                 	19                 	popD_ind5                	popB_ind10               	0.20637398862038292      
DIST(0)          	581          	34                 	20                 	popD_ind5                	popC_ind1                	0.17028612775839083      
DIST(0)          	582          	34                 	21                 	popD_ind5                	popC_ind2                	0.13632496110792791      
DIST(0)          	583          	34                 	22                 	popD_ind5                	popC_ind3                	0.12864368350949146      
DIST(0)          	584          	34                 	23                 	popD_ind5                	popC_ind4                	0.10231511302445452      
DIST(0)          	585          	34                 	24                 	popD_ind5                	popC_ind5                	0.10111958310023154      
DIST(0)          	586          	34                 	25                 	popD_ind5                	popC_ind6                	0.10671319779168148      
DIST(0)          	587          	34                 	26                 	popD_ind5                	popC_ind7                	0.098049834520655618     
DIST(0)          	588          	34                 	27                 	popD_ind5                	popC_ind8                	0.16494759309627052      
DIST(0)          	589          	34                 	28                 	popD_ind5                	popC_ind9                	0.12493239032732519      
DIST(0)          	590          	34                 	29                 	popD_ind5                	popC_ind10               	0.191086241170046        
DIST(0)          	591          	34                 	30                 	popD_ind5                	popD_ind1                	0.10955855872080909      
DIST(0)          	592          	34                 	31                 	popD_ind5                	popD_ind2                	0.063762947474835976     
DIST(0)          	593          	34                 	32                 	popD_ind5                	popD_ind3                	0.09797383948466136      
DIST(0)          	594          	34                 	33                 	popD_ind5                	popD_ind4                	0.070864105330816121     
DIST(0)          	595          	35                 	0                  	popD_ind6                	popA_ind1                	0.18906861575224951      
DIST(0)          	596          	35                 	1                  	popD_ind6                	popA_ind2                	0.17080634713491757      
DIST(0)          	597          	35                 	2                  	popD_ind6                	popA_ind3                	0.15972852315450684      
DIST(0)          	598          	35                 	3                  	popD_ind6                	popA_ind4                	0.19317388248819742      
DIST(0)          	599          	35                 	4                  	popD_ind6                	popA_ind5                	0.21104250591067064      
DIST(0)          	600          	35                 	5                  	popD_ind6                	popA_ind6                	0.23324782242385489      
DIST(0)          	601          	35                 	6                  	popD_ind6                	popA_ind7                	0.21248928031724806      
DIST(0)          	602          	35                 	7                  	popD_ind6                	popA_ind8                	0.19155181086985454      
DIST(0)          	603          	35                 	8                  	popD_ind6                	popA_ind9                	0.17139140749281642      
DIST(0)          	604          	35                 	9                  	popD_ind6                	popA_ind10               	0.1930362464367687       
DIST(0)          	605          	35                 	10                 	popD_ind6                	popB_ind1                	0.1774566876692821       
DIST(0)          	606          	35                 	11                 	popD_ind6                	popB_ind2                	0.16902646191598486      
DIST(0)          	607          	35                 	12                 	popD_ind6                	popB_ind3                	0.22158322180145557      
DIST(0)          	608          	35                 	13                 	popD_ind6                	popB_ind4                	0.20521023961175633      
DIST(0)          	609          	35                 	14                 	popD_ind6                	popB_ind5                	0.24950114791249178      
DIST(0)          	610          	35                 	15                 	popD_ind6                	popB_ind6                	0.27602586879920377      
DIST(0)          	611          	35                 	16                 	popD_ind6                	popB_ind7                	0.17684917109568887      
DIST(0)          	612          	35                 	17                 	popD_ind6                	popB_ind8                	0.1779970924986691       
DIST(0)          	613          	35                 	18                 	popD_ind6                	popB_ind9                	0.21912849783013136      
DIST(0)          	614          	35                 	19                 	popD_ind6                	popB_ind10               	0.19227054846838751      
DIST(0)          	615          	35                 	20                 	popD_ind6                	popC_ind1                	0.16985087038204144      
DIST(0)          	616          	35                 	21                 	popD_ind6                	popC_ind2                	0.13574686785218676      
DIST(0)          	617          	35                 	22                 	popD_ind6                	popC_ind3                	0.14363522563370681      
DIST(0)          	618          	35                 	23                 	popD_ind6                	popC_ind4                	0.14073379222291157      
DIST(0)          	619          	35                 	24                 	popD_ind6                	popC_ind5                	0.15605931194420614      
DIST(0)          	620          	35                 	25                 	popD_ind6                	popC_ind6                	0.13155786527481303      
DIST(0)          	621          	35                 	26                 	popD_ind6                	popC_ind7                	0.1376992537806139       
DIST(0)          	622          	35                 	27                 	popD_ind6                	popC_ind8                	0.17670211034180577      
DIST(0)          	623          	35                 	28                 	popD_ind6                	popC_ind9                	0.15302711426656265      
DIST(0)          	624          	35                 	29                 	popD_ind6                	popC_ind10               	0.20836779159024726      
DIST(0)          	625          	35                 	30                 	popD_ind6                	popD_ind1                	0.058841820962814283     
DIST(0)          	626          	35                 	31                 	popD_ind6                	popD_ind2                	0.085924998986668708     
DIST(0)          	627          	35                 	32                 	popD_ind6                	popD_ind3                	0.05591399061450314      
DIST(0)          	628          	35                 	33                 	popD_ind6                	popD_ind4                	0.054278446686643812     
DIST(0)          	629          	35                 	34                 	popD_ind6                	popD_ind5                	0.083419542736598004     
DIST(0)          	630          	36                 	0                  	popD_ind7                	popA_ind1                	0.21748781832871267      
DIST(0)          	631          	36                 	1                  	popD_ind7                	popA_ind2                	0.18133057044880249      
DIST(0)          	632          	36                 	2                  	popD_ind7                	popA_ind3                	0.15182524495064081      
DIST(0)          	633          	36                 	3                  	popD_ind7                	popA_ind4                	0.22336096806823807      
DIST(0)          	634          	36                 	4                  	popD_ind7                	popA_ind5                	0.23265503843116511      
DIST(0)          	635          	36                 	5                  	popD_ind7                	popA_ind6                	0.20485777633297475      
DIST(0)          	636          	36                 	6                  	popD_ind7                	popA_ind7                	0.24349038518782767      
DIST(0)          	637          	36                 	7                  	popD_ind7                	popA_ind8                	0.18296823283842917      
DIST(0)          	638          	36                 	8                  	popD_ind7                	popA_ind9                	0.15268615997602991      
DIST(0)          	639          	36                 	9                  	popD_ind7                	popA_ind10               	0.2414335536970677       
DIST(0)          	640          	36                 	10                 	popD_ind7                	popB_ind1                	0.23101930380708127      
DIST(0)          	641          	36                 	11                 	popD_ind7                	popB_ind2                	0.12728144468959005      
DIST(0)          	642          	36                 	12                 	popD_ind7                	popB_ind3                	0.15930669423096988      
DIST(0)          	643          	36                 	13                 	popD_ind7                	popB_ind4                	0.22072638825434238      
DIST(0)          	644          	36                 	14                 	popD_ind7                	popB_ind5                	0.25550578838665888      
DIST(0)          	645          	36                 	15                 	popD_ind7                	popB_ind6                	0.2078572933897696       
DIST(0)          	646          	36                 	16                 	popD_ind7                	popB_ind7                	0.2375744362027658       
DIST(0)          	647          	36                 	17                 	popD_ind7                	popB_ind8                	0.21913313280021374      
DIST(0)          	648          	36                 	18                 	popD_ind7                	popB_ind9                	0.24753815852148203      
DIST(0)          	649          	36                 	19                 	popD_ind7                	popB_ind10               	0.2097583578625361       
DIST(0)          	650          	36                 	20                 	popD_ind7                	popC_ind1                	0.16122833726486929      
DIST(0)          	651          	36                 	21                 	popD_ind7                	popC_ind2                	0.13455722684136431      
DIST(0)          	652          	36                 	22                 	popD_ind7                	popC_ind3                	0.15821334044127189      
DIST(0)          	653          	36                 	23                 	popD_ind7                	popC_ind4                	0.17025141661988216      
DIST(0)          	654          	36                 	24                 	popD_ind7                	popC_ind5                	0.13967509150158613      
DIST(0)          	655          	36                 	25                 	popD_ind7                	popC_ind6                	0.17441079092696873      
DIST(0)          	656          	36                 	26                 	popD_ind7                	popC_ind7                	0.16216400083428656      
DIST(0)          	657          	36                 	27                 	popD_ind7                	popC_ind8                	0.17065790589462435      
DIST(0)          	658          	36                 	28                 	popD_ind7                	popC_ind9                	0.19962883261429853      
DIST(0)          	659          	36                 	29                 	popD_ind7                	popC_ind10               	0.18540859914235105      
DIST(0)          	660          	36                 	30                 	popD_ind7                	popD_ind1                	0.081386439574588465     
DIST(0)          	661          	36                 	31                 	popD_ind7                	popD_ind2                	0.059289016191081186     
DIST(0)          	662          	36                 	32                 	popD_ind7                	popD_ind3                	0.066235445652526495     
DIST(0)          	663          	36                 	33                 	popD_ind7                	popD_ind4                	0.052475010496292815     
DIST(0)          	664          	36                 	34                 	popD_ind7                	popD_ind5                	0.10335837115998309      
DIST(0)          	665          	36                 	35                 	popD_ind7                	popD_ind6                	0.10181903563687017      
DIST(0)          	666          	37                 	0                  	popD_ind8                	popA_ind1                	0.28611408867877519      
DIST(0)          	667          	37                 	1                  	popD_ind8                	popA_ind2                	0.19701157640374312      
DIST(0)          	668          	37                 	2                  	popD_ind8                	popA_ind3                	0.20851366427915063      
DIST(0)          	669          	37                 	3                  	popD_ind8                	popA_ind4                	0.2175316212354968       
DIST(0)          	670          	37                 	4                  	popD_ind8                	popA_ind5                	0.23075670961912637      
DIST(0)          	671          	37                 	5                  	popD_ind8                	popA_ind6                	0.26979832917824487      
DIST(0)          	672          	37                 	6                  	popD_ind8                	popA_ind7                	0.24279219748868724      
DIST(0)          	673          	37                 	7                  	popD_ind8                	popA_ind8                	0.18677507061986889      
DIST(0)          	674          	37                 	8                  	popD_ind8                	popA_ind9                	0.23162921250311166      
DIST(0)          	675          	37                 	9                  	popD_ind8                	popA_ind10               	0.26385518040665779      
DIST(0)          	676          	37                 	10                 	popD_ind8                	popB_ind1                	0.23180252526263262      
DIST(0)          	677          	37                 	11                 	popD_ind8                	popB_ind2                	0.2564651466160206       
DIST(0)          	678          	37                 	12                 	popD_ind8                	popB_ind3                	0.22949421606864762      
DIST(0)          	679          	37                 	13                 	popD_ind8                	popB_ind4                	0.21021289975470153      
DIST(0)          	680          	37                 	14                 	popD_ind8                	popB_ind5                	0.24035674846257954      
DIST(0)          	681          	37                 	15                 	popD_ind8                	popB_ind6                	0.23246961300340269      
DIST(0)          	682          	37                 	16                 	popD_ind8                	popB_ind7                	0.24919367820131316      
DIST(0)          	683          	37                 	17                 	popD_ind8                	popB_ind8                	0.23188610409305635      
DIST(0)          	684          	37                 	18                 	popD_ind8                	popB_ind9                	0.21870880054468211      
DIST(0)          	685          	37                 	19                 	popD_ind8                	popB_ind10               	0.2750491602878905       
DIST(0)          	686          	37                 	20                 	popD_ind8                	popC_ind1                	0.14123894476618121      
DIST(0)          	687          	37                 	21                 	popD_ind8                	popC_ind2                	0.14543631086115855      
DIST(0)          	688          	37                 	22                 	popD_ind8                	popC_ind3                	0.14769445662723027      
DIST(0)          	689          	37                 	23                 	popD_ind8                	popC_ind4                	0.1324365005803228       
DIST(0)          	690          	37                 	24                 	popD_ind8                	popC_ind5                	0.13959193467527092      
DIST(0)          	691          	37                 	25                 	popD_ind8                	popC_ind6                	0.18151937177781793      
DIST(0)          	692          	37                 	26                 	popD_ind8                	popC_ind7                	0.21248647886909616      
DIST(0)          	693          	37                 	27                 	popD_ind8                	popC_ind8                	0.13484168233253779      
DIST(0)          	694          	37                 	28                 	popD_ind8                	popC_ind9                	0.15687661938295838      
DIST(0)          	695          	37                 	29                 	popD_ind8                	popC_ind10               	0.12457208959514103      
DIST(0)          	696          	37                 	30                 	popD_ind8                	popD_ind1                	0.078822163626686204     
DIST(0)          	697          	37                 	31                 	popD_ind8                	popD_ind2                	0.031588077619131037     
DIST(0)          	698          	37                 	32                 	popD_ind8                	popD_ind3                	0.11223839620790221      
DIST(0)          	699          	37                 	33                 	popD_ind8                	popD_ind4                	0.10120119662321263      
DIST(0)          	700          	37                 	34                 	popD_ind8                	popD_ind5                	0.11154495870553233      
DIST(0)          	701          	37                 	35                 	popD_ind8                	popD_ind6                	0.12718764441321856      
DIST(0)          	702          	37                 	36                 	popD_ind8                	popD_ind7                	0.080805289358370724     
DIST(0)          	703          	38                 	0                  	popD_ind9                	popA_ind1                	0.19386377200036506      
DIST(0)          	704          	38                 	1                  	popD_ind9                	popA_ind2                	0.17041531683237982      
DIST(0)          	705          	38                 	2                  	popD_ind9                	popA_ind3                	0.1801110486003773       
DIST(0)          	706          	38                 	3                  	popD_ind9                	popA_ind4                	0.18727158297348423      
DIST(0)          	707          	38                 	4                  	popD_ind9                	popA_ind5                	0.17363076858038598      
DIST(0)          	708          	38                 	5                  	popD_ind9                	popA_ind6                	0.18137935721932291      
DIST(0)          	709          	38                 	6                  	popD_ind9                	popA_ind7                	0.15964983626200088      
DIST(0)          	710          	38                 	7                  	popD_ind9                	popA_ind8                	0.19857084662885982      
DIST(0)          	711          	38                 	8                  	popD_ind9                	popA_ind9                	0.1896546613283232       
DIST(0)          	712          	38                 	9                  	popD_ind9                	popA_ind10               	0.1868457516815038       
DIST(0)          	713          	38                 	10                 	popD_ind9                	popB_ind1                	0.21999460631057782      
DIST(0)          	714          	38                 	11                 	popD_ind9                	popB_ind2                	0.160897358214598        
DIST(0)          	715          	38                 	12                 	popD_ind9                	popB_ind3                	0.16183356523384151      
DIST(0)          	716          	38                 	13                 	popD_ind9                	popB_ind4                	0.17926301699703356      
DIST(0)          	717          	38                 	14                 	popD_ind9                	popB_ind5                	0.20773487075049935      
DIST(0)          	718          	38                 	15                 	popD_ind9                	popB_ind6                	0.16792936705061307      
DIST(0)          	719          	38                 	16                 	popD_ind9                	popB_ind7                	0.1935185309877249       
DIST(0)          	720          	38                 	17                 	popD_ind9                	popB_ind8                	0.19781970164311488      
DIST(0)          	721          	38                 	18                 	popD_ind9                	popB_ind9                	0.20858438752073355      
DIST(0)          	722          	38                 	19                 	popD_ind9                	popB_ind10               	0.24375939834005073      
DIST(0)          	723          	38                 	20                 	popD_ind9                	popC_ind1                	0.14987999185384393      
DIST(0)          	724          	38                 	21                 	popD_ind9                	popC_ind2                	0.1090137443288704       
DIST(0)          	725          	38                 	22                 	popD_ind9                	popC_ind3                	0.17818709553400197      
DIST(0)          	726          	38                 	23                 	popD_ind9                	popC_ind4                	0.1145785571156178       
DIST(0)          	727          	38                 	24                 	popD_ind9                	popC_ind5                	0.10746869968649356      
DIST(0)          	728          	38                 	25                 	popD_ind9                	popC_ind6                	0.13536378475492475      
DIST(0)          	729          	38                 	26                 	popD_ind9                	popC_ind7                	0.19869567622299122      
DIST(0)          	730          	38                 	27                 	popD_ind9                	popC_ind8                	0.12890692492276618      
DIST(0)          	731          	38                 	28                 	popD_ind9                	popC_ind9                	0.10481242462390895      
DIST(0)          	732          	38                 	29                 	popD_ind9                	popC_ind10               	0.11043098197863575      
DIST(0)          	733          	38                 	30                 	popD_ind9                	popD_ind1                	0.085772701730966197     
DIST(0)          	734          	38                 	31                 	popD_ind9                	popD_ind2                	0.061212496725935248     
DIST(0)          	735          	38                 	32                 	popD_ind9                	popD_ind3                	0.1008486333597248       
DIST(0)          	736          	38                 	33                 	popD_ind9                	popD_ind4                	0.044549299427491008     
DIST(0)          	737          	38                 	34                 	popD_ind9                	popD_ind5                	0.13737949274062231      
DIST(0)          	738          	38                 	35                 	popD_ind9                	popD_ind6                	0.11763513014769206      
DIST(0)          	739          	38                 	36                 	popD_ind9                	popD_ind7                	0.073765825390509548     
DIST(0)          	740          	38                 	37                 	popD_ind9                	popD_ind8                	0.058150291609240713     
DIST(0)          	741          	39                 	0                  	popD_ind10               	popA_ind1                	0.23591998461270308      
DIST(0)          	742          	39                 	1                  	popD_ind10               	popA_ind2                	0.18141088869471153      
DIST(0)          	743          	39                 	2                  	popD_ind10               	popA_ind3                	0.1751242271860366       
DIST(0)          	744          	39                 	3                  	popD_ind10               	popA_ind4                	0.21880728377739106      
DIST(0)          	745          	39                 	4                  	popD_ind10               	popA_ind5                	0.23389374940547453      
DIST(0)          	746          	39                 	5                  	popD_ind10               	popA_ind6                	0.25865775635084853      
DIST(0)          	747          	39                 	6                  	popD_ind10               	popA_ind7                	0.25736482074151124      
DIST(0)          	748          	39                 	7                  	popD_ind10               	popA_ind8                	0.22586486844133075      
DIST(0)          	749          	39                 	8                  	popD_ind10               	popA_ind9                	0.1981339353129731       
DIST(0)          	750          	39                 	9                  	popD_ind10               	popA_ind10               	0.24799146793248905      
DIST(0)          	751          	39                 	10                 	popD_ind10               	popB_ind1                	0.26950201588950762      
DIST(0)          	752          	39                 	11                 	popD_ind10               	popB_ind2                	0.18603445591209461      
DIST(0)          	753          	39                 	12                 	popD_ind10               	popB_ind3                	0.22598027326955611      
DIST(0)          	754          	39                 	13                 	popD_ind10               	popB_ind4                	0.2442043550545287       
DIST(0)          	755          	39                 	14                 	popD_ind10               	popB_ind5                	0.26915977775154032      
DIST(0)          	756          	39                 	15                 	popD_ind10               	popB_ind6                	0.26054458549968323      
DIST(0)          	757          	39                 	16                 	popD_ind10               	popB_ind7                	0.23159659610687203      
DIST(0)          	758          	39                 	17                 	popD_ind10               	popB_ind8                	0.25040070451259233      
DIST(0)          	759          	39                 	18                 	popD_ind10               	popB_ind9                	0.22337389092966381      
DIST(0)          	760          	39                 	19                 	popD_ind10               	popB_ind10               	0.2053796552416069       
DIST(0)          	761          	39                 	20                 	popD_ind10               	popC_ind1                	0.18979746442769113      
DIST(0)          	762          	39                 	21                 	popD_ind10               	popC_ind2                	0.14327432104873153      
DIST(0)          	763          	39                 	22                 	popD_ind10               	popC_ind3                	0.19537235831398664      
DIST(0)          	764          	39                 	23                 	popD_ind10               	popC_ind4                	0.14155429190999963      
DIST(0)          	765          	39                 	24                 	popD_ind10               	popC_ind5                	0.12936714956764023      
DIST(0)          	766          	39                 	25                 	popD_ind10               	popC_ind6                	0.13009028477779364      
DIST(0)          	767          	39                 	26                 	popD_ind10               	popC_ind7                	0.15849379761903853      
DIST(0)          	768          	39                 	27                 	popD_ind10               	popC_ind8                	0.13629876919713782      
DIST(0)          	769          	39                 	28                 	popD_ind10               	popC_ind9                	0.15005500084582213      
DIST(0)          	770          	39                 	29                 	popD_ind10               	popC_ind10               	0.19704344649624517      
DIST(0)          	771          	39                 	30                 	popD_ind10               	popD_ind1                	0.095958891846744498     
DIST(0)          	772          	39                 	31                 	popD_ind10               	popD_ind2                	0.076421736302500723     
DIST(0)          	773          	39                 	32                 	popD_ind10               	popD_ind3                	0.076350578502792574     
DIST(0)          	774          	39                 	33                 	popD_ind10               	popD_ind4                	0.10395933275132016      
DIST(0)          	775          	39                 	34                 	popD_ind10               	popD_ind5                	0.075444102704649463     
DIST(0)          	776          	39                 	35                 	popD_ind10               	popD_ind6                	0.097702750023059115     
DIST(0)          	777          	39                 	36                 	popD_ind10               	popD_ind7                	0.09300147959056855      
DIST(0)          	778          	39                 	37                 	popD_ind10               	popD_ind8                	0.070368147500554382     
DIST(0)          	779          	39                 	38                 	popD_ind10               	popD_ind9                	0.1273092014238813       
# DMAT(MATRIX_INDEX), Distance matrix in matrix format:
# DMAT(0)        	[2]ITEMS                 	popA_ind1                	popA_ind2                	popA_ind3                	popA_ind4                	popA_ind5                	popA_ind6                	popA_ind7                	popA_ind8                	popA_ind9                	popA_ind10               	popB_ind1                	popB_ind2                	popB_ind3                	popB_ind4                	popB_ind5                	popB_ind6                	popB_ind7                	popB_ind8                	popB_ind9                	popB_ind10               	popC_ind1                	popC_ind2                	popC_ind3                	popC_ind4                	popC_ind5                	popC_ind6                	popC_ind7                	popC_ind8                	popC_ind9                	popC_ind10               	popD_ind1                	popD_ind2                	popD_ind3                	popD_ind4                	popD_ind5                	popD_ind6                	popD_ind7                	popD_ind8                	popD_ind9                	popD_ind10               
DMAT(0)          	popA_ind1                	
DMAT(0)          	popA_ind2                	0.11467617878074979      
DMAT(0)          	popA_ind3                	0.093092202947423874     	0.043429845634977288     
DMAT(0)          	popA_ind4                	0.09906854300806188      	0.14476271136798249      	0.09250928862574137      
DMAT(0)          	popA_ind5                	0.11548749475140634      	0.10235826812419881      	0.095219871094828198     	0.10924768687545443      
DMAT(0)          	popA_ind6                	0.10717339691795179      	0.081245857649091494     	0.10665280767883598      	0.0992253751337005       	0.084457781887820244     
DMAT(0)          	popA_ind7                	0.10530027483526966      	0.10952915292971824      	0.056029316909226953     	0.11165936466273685      	0.11752068019191365      	0.099017757040388107     
DMAT(0)          	popA_ind8                	0.097166560252935227     	0.087725183305922882     	0.035463864283368375     	0.10971550268283765      	0.052211145413247576     	0.080090237697028069     	0.081349577339648454     
DMAT(0)          	popA_ind9                	0.090441910531874212     	0.095867900816094878     	0.058549180362201261     	0.1026648719420762       	0.069954678492391426     	0.092754438818220819     	0.067255283583122355     	0.081902410789416374     
DMAT(0)          	popA_ind10               	0.12644992674147426      	0.096636012476263314     	0.10716580364420708      	0.10184357392425492      	0.10530742687306804      	0.079663031092544523     	0.13024246593697633      	0.08651233086631463      	0.091715136836088701     
DMAT(0)          	popB_ind1                	0.2336180839991483       	0.1895475907764064       	0.2028394442494629       	0.23700842812959011      	0.21720537528782002      	0.16646532629003258      	0.1794056680155674       	0.18782183744554481      	0.17880978453375057      	0.13315769854599507      
DMAT(0)          	popB_ind2                	0.15574930521351982      	0.12902711380341364      	0.18546996912365243      	0.2009515144952744       	0.13921199302111537      	0.1186383060851682       	0.2186548244483244       	0.16953188761613003      	0.13355896813518409      	0.13685061378832791      	0.082160885097566222     
DMAT(0)          	popB_ind3                	0.18435734223985237      	0.18073172098884421      	0.12921894115135346      	0.12284745856950652      	0.18104606093364706      	0.19471725344017071      	0.19364555475485268      	0.11517994569806497      	0.1348604073845501       	0.19398220480765638      	0.16274535265075507      	0.13373797837552712      
DMAT(0)          	popB_ind4                	0.13365602036443922      	0.17171434159611629      	0.14385416636385348      	0.13620392951227267      	0.13844854551032743      	0.20507900555218073      	0.16674800759215752      	0.19158365508898345      	0.14481723658856765      	0.17216895148979988      	0.10146661607076186      	0.12618862661658062      	0.058739516198463057     
DMAT(0)          	popB_ind5                	0.21820072728406747      	0.1465067556144577       	0.19196571222817477      	0.1761893161168778       	0.16128517489385935      	0.15657767928982227      	0.17567679918954104      	0.16873546816436635      	0.19334619577350509      	0.14622659497598675      	0.076306207154471978     	0.066529877811036531     	0.1409409842754705       	0.096059387487813133     
DMAT(0)          	popB_ind6                	0.13848934374173905      	0.15612309594158841      	0.17415526086063793      	0.19717244117112806      	0.20189920751847495      	0.21586937475895537      	0.18762449822932767      	0.19803684293859386      	0.17157542180045349      	0.152054208614887        	0.1370090321125082       	0.087766385293120752     	0.14859284705566383      	0.13156715752246437      	0.10742285964688827      
DMAT(0)          	popB_ind7                	0.17973240429294551      	0.14078806644125302      	0.091761684712605401     	0.14915398801606838      	0.16024652644638926      	0.16262642621955026      	0.12382788967982038      	0.12664600605320159      	0.11738137106945715      	0.15210423771205939      	0.082060786862680438     	0.13608893161450786      	0.080284930300773313     	0.069619761483182135     	0.13062814343793533      	0.14559357684436564      
DMAT(0)          	popB_ind8                	0.16567511205365193      	0.12640773566882468      	0.1708866189778657       	0.23133671988418389      	0.19225338752473922      	0.19726913409919919      	0.1479487393021007       	0.17227982719132304      	0.14959957714123032      	0.15249983598143485      	0.096185128202881956     	0.12356963919485955      	0.097396995678103115     	0.1170702575416423       	0.077202720953757881     	0.12856793997731342      	0.12011526516704371      
DMAT(0)          	popB_ind9                	0.12538366546380877      	0.1536787031523098       	0.10515210094250803      	0.15781540012958445      	0.12211001480676476      	0.13742681971106568      	0.11147969497017542      	0.14821969503273411      	0.095865949319817401     	0.15796362529584287      	0.11159237196376758      	0.10566166338278885      	0.1183479934870411       	0.097778578496353533     	0.092115999449317809     	0.11123602906816318      	0.079927548958855121     	0.097300091716630593     
DMAT(0)          	popB_ind10               	0.2130766619311581       	0.17037706765265032      	0.14639700483239165      	0.19214745347349027      	0.13061699069084559      	0.15657229526571761      	0.1698556281188352       	0.18171893336639938      	0.16719664771957576      	0.1634737129064138       	0.12687343944407434      	0.1443034076029559       	0.071379753531828224     	0.074913241182063398     	0.11966181426920425      	0.14635490020317696      	0.093136998894170486     	0.079366584664830098     	0.1142753889662504       
DMAT(0)          	popC_ind1                	0.2015332824908172       	0.21319717169203276      	0.1993509842944885       	0.22274743747405046      	0.13482631394144037      	0.22349667533920586      	0.21226372044518818      	0.15416149791707109      	0.14624866514036186      	0.19273603119366231      	0.18992584418755429      	0.13916576111370943      	0.218718029433946        	0.23894372190687277      	0.17615682862512977      	0.21229738138034671      	0.25276437247717165      	0.22832392155407685      	0.20775044417262856      	0.1441944134983002       
DMAT(0)          	popC_ind2                	0.22010928878735495      	0.19596048314884856      	0.13348044864677658      	0.22768146183832755      	0.17210693161572244      	0.22974283767662557      	0.20868695076057758      	0.21980391670649729      	0.20436847500812053      	0.23773925327304762      	0.2411019389417306       	0.18274924278785148      	0.22482126011468145      	0.21303988719470016      	0.24704618511131388      	0.17241113608750086      	0.22743761668678628      	0.24500901223963262      	0.17401907242445552      	0.17641728721519201      	0.085929952914894553     
DMAT(0)          	popC_ind3                	0.20992072881488485      	0.19856025357945112      	0.18137186054747748      	0.22563648648046836      	0.19355061640016907      	0.20568966173104586      	0.18097502863792983      	0.19008316661633079      	0.20914506792198181      	0.196616965692434        	0.20516226119894759      	0.17419158370645471      	0.22366785596444297      	0.23600090782376576      	0.24198609467351978      	0.19888742525337738      	0.21567635282876946      	0.20019062125200141      	0.20930096488445013      	0.19961998820517368      	0.064587240380522284     	0.083556136378852702     
DMAT(0)          	popC_ind4                	0.18011973954276256      	0.17799342432036944      	0.14461749966416013      	0.19707754359541055      	0.081157311261563433     	0.19920225108790207      	0.23824256705359415      	0.12737489517049935      	0.13183984482082817      	0.16902390438926701      	0.16848696867103974      	0.14855713001276011      	0.22618623471090815      	0.19733227421567026      	0.219284944468256        	0.14381946282788433      	0.18400919938333199      	0.21535393907604786      	0.17428643282842965      	0.16779312578307787      	0.088996394857584024     	0.063669140840762459     	0.087365595598285273     
DMAT(0)          	popC_ind5                	0.20325323130980855      	0.14702450642836398      	0.12378919856248134      	0.16490775905017863      	0.12756045244004605      	0.1527573957747721       	0.14767875278433937      	0.11134245564017317      	0.15457974651569453      	0.16103514198867447      	0.18153955539439673      	0.12645621735530055      	0.21219604150214655      	0.22434630118986329      	0.21254986108129215      	0.16881940483549152      	0.14723056087689204      	0.20735358341427793      	0.15812908243462714      	0.12767489028796095      	0.038373039060623441     	0.092059644589807504     	0.074413818107818733     	0.053726777929373533     
DMAT(0)          	popC_ind6                	0.20984031673621631      	0.17442468605602207      	0.2589571139309198       	0.2116313684226776       	0.16355731741209178      	0.26326882256131295      	0.26393090198139946      	0.16815725826357403      	0.22714663644873606      	0.24610165008253948      	0.24285448097296569      	0.18794323710484434      	0.24428109312872376      	0.18087135708310498      	0.18653439343397094      	0.20582109677591914      	0.22320203943759706      	0.21083984955048846      	0.18504443247286206      	0.1909788577390174       	0.061336527554696936     	0.060911166662707003     	0.088665456548880486     	0.074081134186570141     	0.053431671951617909     
DMAT(0)          	popC_ind7                	0.23655857744765843      	0.24830740085268679      	0.2216585993132843       	0.18967820697965593      	0.2013899855208609       	0.28795560981523416      	0.26873541148604974      	0.19675588091388391      	0.22811461996768145      	0.20867271024605338      	0.24061109001033371      	0.19877923217784799      	0.21547008840665582      	0.19767769172169311      	0.27509638561703642      	0.25374463338126252      	0.17857892695742406      	0.22978234231581426      	0.22425103724399484      	0.19836430073578473      	0.080468945192271019     	0.077573278036358176     	0.096624243274473623     	0.078420141511245114     	0.085692851086862684     	0.094015717217159234     
DMAT(0)          	popC_ind8                	0.23171725978169411      	0.22090372507239436      	0.20809777426540216      	0.19897184497235865      	0.2098146614756794       	0.23941607366271042      	0.20018477116636771      	0.17347278125565394      	0.20419846762542376      	0.21781060486206261      	0.22480874276002177      	0.27709420638695614      	0.22787476566118275      	0.19083714046180256      	0.30645922919112667      	0.23059273142635181      	0.23069574344760835      	0.26584523918957537      	0.22086612680901901      	0.19970167979323389      	0.038367044363519641     	0.057603914036350262     	0.10002365391216407      	0.074811407190030099     	0.049541727827786639     	0.083376674419933378     	0.10231720437086898      
DMAT(0)          	popC_ind9                	0.27893889901517865      	0.21410309795138521      	0.19286296697013489      	0.18686749308825451      	0.15456903913146047      	0.21061187153154842      	0.22837868173009898      	0.1637846029129669       	0.17143383683045679      	0.14333689465982716      	0.2403157110690059       	0.19753342149877423      	0.24909517760238462      	0.21178624032341026      	0.23701875375463244      	0.20842888691417302      	0.27955586354207185      	0.29242445900054809      	0.25257668494528768      	0.19507042811246808      	0.062541815271601547     	0.096198048711564468     	0.1280191774205563       	0.095481720882076723     	0.046700724870551005     	0.078023295539278489     	0.094692791475228327     	0.094952695531240541     
DMAT(0)          	popC_ind10               	0.26359692965102521      	0.24990291569905126      	0.24506414677891708      	0.27575413391200859      	0.24426070118934098      	0.22138250336297252      	0.22809973495139552      	0.22615354130733101      	0.24653119815044355      	0.2294934646127805       	0.20338124312427919      	0.23936921589388946      	0.24991241331889907      	0.22499187380586272      	0.27313820987458215      	0.25544194022353872      	0.2754649220079905       	0.2684169366788503       	0.21351601714660881      	0.19710585973442807      	0.084133763987661189     	0.089002603093127908     	0.090197021640572694     	0.070476984843124091     	0.069405838498509606     	0.081919128661087498     	0.10398249075452988      	0.088913855721502233     	0.061781873583093873     
DMAT(0)          	popD_ind1                	0.1557687437607343       	0.18354603047325468      	0.15548244825433166      	0.23636189166322444      	0.16065407838892959      	0.25146596833430723      	0.18426987979109452      	0.1424256545264975       	0.18239425024363842      	0.20453525675933962      	0.20477225865304255      	0.20888693041062814      	0.19351193017586943      	0.199913543115711        	0.24126918206266978      	0.17394113001215428      	0.18805765501821292      	0.18541141583156168      	0.15971641445464083      	0.20158372573305072      	0.12984999798110172      	0.093924844356783055     	0.14958527171031463      	0.087817162380279215     	0.074551563154910017     	0.11993332403573877      	0.16517261220996154      	0.11194600643558583      	0.13796409469097348      	0.16068347267613489      
DMAT(0)          	popD_ind2                	0.2252919853265799       	0.15273144787650705      	0.23137037387074691      	0.21336392171997742      	0.18413506961360532      	0.20093618363776045      	0.17424317963019315      	0.22066812525444179      	0.16343172198789713      	0.2683787525760486       	0.2176982159169541       	0.22592371461100735      	0.27272862947566057      	0.20245036636524577      	0.25020944057727357      	0.26265371451431879      	0.2265333001434483       	0.23948524980857971      	0.21237658573738638      	0.19899818404453437      	0.16846483892402273      	0.166728731623842        	0.16175438219767974      	0.082048556450437959     	0.065738464466928354     	0.12302478124949576      	0.1624192388846194       	0.17524499344050826      	0.11883894416167089      	0.17443339134304406      	0.064576161974153279     
DMAT(0)          	popD_ind3                	0.22101561676569037      	0.091315068731080598     	0.17047831749440684      	0.18281749583573956      	0.19268107311000748      	0.22939261178266604      	0.17338542867880857      	0.19837001276276728      	0.15008270114662817      	0.18563201568799226      	0.21893510830319929      	0.18349102271350159      	0.21277993095707076      	0.17493860817319762      	0.19796806509855142      	0.17308135954369908      	0.20618120362975073      	0.1463432015373973       	0.17480534081730839      	0.17374194623019329      	0.10930288882255604      	0.10846395089436804      	0.12783401100466404      	0.11469374835607198      	0.081969681336522646     	0.13636551983836287      	0.16349176700277329      	0.11927721961442059      	0.12640741666249491      	0.20027109942585425      	0.055549174616355043     	0.084741591710787309     
DMAT(0)          	popD_ind4                	0.20867651590808914      	0.1586076732132409       	0.20642891361451218      	0.19888747055935219      	0.18392410152495564      	0.22144594242960736      	0.17902726726923524      	0.20990936913047709      	0.19507804468020856      	0.18038684264531946      	0.23919289155923101      	0.15478947561292911      	0.19941233218382945      	0.20573764340109796      	0.25391447973445547      	0.23530740789842947      	0.22514748421643299      	0.21471598894825331      	0.19174464449363682      	0.21651897983173204      	0.14684020519339983      	0.12154685529474632      	0.16840734321912113      	0.1310283516519139       	0.10528296756241991      	0.092455188594602639     	0.15531964686520139      	0.17810137419091548      	0.11642410183818176      	0.16329440387638294      	0.040291159308572445     	0.071182101174123089     	0.050850306031042618     
DMAT(0)          	popD_ind5                	0.2288501073715965       	0.23106126738994098      	0.21736566698621612      	0.21130195979700117      	0.25371929719601982      	0.25158744905931074      	0.23206171111174118      	0.22433396965282579      	0.20455029119506379      	0.22525837630830783      	0.22001406799952861      	0.17962672021236251      	0.19070710283874676      	0.17960029739689223      	0.24667554847277942      	0.24720539059618932      	0.22253210946423452      	0.22024858207637088      	0.22087491212067947      	0.20637398862038292      	0.17028612775839083      	0.13632496110792791      	0.12864368350949146      	0.10231511302445452      	0.10111958310023154      	0.10671319779168148      	0.098049834520655618     	0.16494759309627052      	0.12493239032732519      	0.191086241170046        	0.10955855872080909      	0.063762947474835976     	0.09797383948466136      	0.070864105330816121     
DMAT(0)          	popD_ind6                	0.18906861575224951      	0.17080634713491757      	0.15972852315450684      	0.19317388248819742      	0.21104250591067064      	0.23324782242385489      	0.21248928031724806      	0.19155181086985454      	0.17139140749281642      	0.1930362464367687       	0.1774566876692821       	0.16902646191598486      	0.22158322180145557      	0.20521023961175633      	0.24950114791249178      	0.27602586879920377      	0.17684917109568887      	0.1779970924986691       	0.21912849783013136      	0.19227054846838751      	0.16985087038204144      	0.13574686785218676      	0.14363522563370681      	0.14073379222291157      	0.15605931194420614      	0.13155786527481303      	0.1376992537806139       	0.17670211034180577      	0.15302711426656265      	0.20836779159024726      	0.058841820962814283     	0.085924998986668708     	0.05591399061450314      	0.054278446686643812     	0.083419542736598004     
DMAT(0)          	popD_ind7                	0.21748781832871267      	0.18133057044880249      	0.15182524495064081      	0.22336096806823807      	0.23265503843116511      	0.20485777633297475      	0.24349038518782767      	0.18296823283842917      	0.15268615997602991      	0.2414335536970677       	0.23101930380708127      	0.12728144468959005      	0.15930669423096988      	0.22072638825434238      	0.25550578838665888      	0.2078572933897696       	0.2375744362027658       	0.21913313280021374      	0.24753815852148203      	0.2097583578625361       	0.16122833726486929      	0.13455722684136431      	0.15821334044127189      	0.17025141661988216      	0.13967509150158613      	0.17441079092696873      	0.16216400083428656      	0.17065790589462435      	0.19962883261429853      	0.18540859914235105      	0.081386439574588465     	0.059289016191081186     	0.066235445652526495     	0.052475010496292815     	0.10335837115998309      	0.10181903563687017      
DMAT(0)          	popD_ind8                	0.28611408867877519      	0.19701157640374312      	0.20851366427915063      	0.2175316212354968       	0.23075670961912637      	0.26979832917824487      	0.24279219748868724      	0.18677507061986889      	0.23162921250311166      	0.26385518040665779      	0.23180252526263262      	0.2564651466160206       	0.22949421606864762      	0.21021289975470153      	0.24035674846257954      	0.23246961300340269      	0.24919367820131316      	0.23188610409305635      	0.21870880054468211      	0.2750491602878905       	0.14123894476618121      	0.14543631086115855      	0.14769445662723027      	0.1324365005803228       	0.13959193467527092      	0.18151937177781793      	0.21248647886909616      	0.13484168233253779      	0.15687661938295838      	0.12457208959514103      	0.078822163626686204     	0.031588077619131037     	0.11223839620790221      	0.10120119662321263      	0.11154495870553233      	0.12718764441321856      	0.080805289358370724     
DMAT(0)          	popD_ind9                	0.19386377200036506      	0.17041531683237982      	0.1801110486003773       	0.18727158297348423      	0.17363076858038598      	0.18137935721932291      	0.15964983626200088      	0.19857084662885982      	0.1896546613283232       	0.1868457516815038       	0.21999460631057782      	0.160897358214598        	0.16183356523384151      	0.17926301699703356      	0.20773487075049935      	0.16792936705061307      	0.1935185309877249       	0.19781970164311488      	0.20858438752073355      	0.24375939834005073      	0.14987999185384393      	0.1090137443288704       	0.17818709553400197      	0.1145785571156178       	0.10746869968649356      	0.13536378475492475      	0.19869567622299122      	0.12890692492276618      	0.10481242462390895      	0.11043098197863575      	0.085772701730966197     	0.061212496725935248     	0.1008486333597248       	0.044549299427491008     	0.13737949274062231      	0.11763513014769206      	0.073765825390509548     	0.058150291609240713     
DMAT(0)          	popD_ind10               	0.23591998461270308      	0.18141088869471153      	0.1751242271860366       	0.21880728377739106      	0.23389374940547453      	0.25865775635084853      	0.25736482074151124      	0.22586486844133075      	0.1981339353129731       	0.24799146793248905      	0.26950201588950762      	0.18603445591209461      	0.22598027326955611      	0.2442043550545287       	0.26915977775154032      	0.26054458549968323      	0.23159659610687203      	0.25040070451259233      	0.22337389092966381      	0.2053796552416069       	0.18979746442769113      	0.14327432104873153      	0.19537235831398664      	0.14155429190999963      	0.12936714956764023      	0.13009028477779364      	0.15849379761903853      	0.13629876919713782      	0.15005500084582213      	0.19704344649624517      	0.095958891846744498     	0.076421736302500723     	0.076350578502792574     	0.10395933275132016      	0.075444102704649463     	0.097702750023059115     	0.09300147959056855      	0.070368147500554382     	0.1273092014238813       
## -> Following information is specific to distance matrix with MATRIX_INDEX=1 in this file:
# DIST(MATRIX_INDEX), Pairwise distances in the distance matrix:
# DIST(1)        	[4]PAIR_INDEX	[5]FIRST_ITEM_INDEX	[6]SECOND_ITEM_INDEX	[7]FIRST_ITEM_NAME       	[8]SECOND_ITEM_NAME      	[9]DISTANCE              
DIST(1)          	0            	1                  	0                  	popA_ind2                	popA_ind1                	0.11467617878074976      
DIST(1)          	1            	2                  	0                  	popA_ind3                	popA_ind1                	0.093092202947423805     
DIST(1)          	2            	2                  	1                  	popA_ind3                	popA_ind2                	0.043429845634977267     
DIST(1)          	3            	3                  	0                  	popA_ind4                	popA_ind1                	0.099068543008061838     
DIST(1)          	4            	3                  	1                  	popA_ind4                	popA_ind2                	0.14476271136798247      
DIST(1)          	5            	3                  	2                  	popA_ind4                	popA_ind3                	0.092509288625741343     
DIST(1)          	6            	4                  	0                  	popA_ind5                	popA_ind1                	0.11548749475140628      
DIST(1)          	7            	4                  	1                  	popA_ind5                	popA_ind2                	0.10235826812419882      
DIST(1)          	8            	4                  	2                  	popA_ind5                	popA_ind3                	0.09521987109482824      
DIST(1)          	9            	4                  	3                  	popA_ind5                	popA_ind4                	0.10924768687545468      
DIST(1)          	10           	5                  	0                  	popA_ind6                	popA_ind1                	0.10717339691795183      
DIST(1)          	11           	5                  	1                  	popA_ind6                	popA_ind2                	0.08124585764909141      
DIST(1)          	12           	5                  	2                  	popA_ind6                	popA_ind3                	0.10665280767883602      
DIST(1)          	13           	5                  	3                  	popA_ind6                	popA_ind4                	0.099225375133700472     
DIST(1)          	14           	5                  	4                  	popA_ind6                	popA_ind5                	0.084457781887820299     
DIST(1)          	15           	6                  	0                  	popA_ind7                	popA_ind1                	0.10530027483526966      
DIST(1)          	16           	6                  	1                  	popA_ind7                	popA_ind2                	0.10952915292971828      
DIST(1)          	17           	6                  	2                  	popA_ind7                	popA_ind3                	0.056029316909226939     
DIST(1)          	18           	6                  	3                  	popA_ind7                	popA_ind4                	0.11165936466273696      
DIST(1)          	19           	6                  	4                  	popA_ind7                	popA_ind5                	0.11752068019191368      
DIST(1)          	20           	6                  	5                  	popA_ind7                	popA_ind6                	0.099017757040388121     
DIST(1)          	21           	7                  	0                  	popA_ind8                	popA_ind1                	0.097166560252935116     
DIST(1)          	22           	7                  	1                  	popA_ind8                	popA_ind2                	0.087725183305923007     
DIST(1)          	23           	7                  	2                  	popA_ind8                	popA_ind3                	0.035463864283368333     
DIST(1)          	24           	7                  	3                  	popA_ind8                	popA_ind4                	0.10971550268283756      
DIST(1)          	25           	7                  	4                  	popA_ind8                	popA_ind5                	0.052211145413247562     
DIST(1)          	26           	7                  	5                  	popA_ind8                	popA_ind6                	0.080090237697028083     
DIST(1)          	27           	7                  	6                  	popA_ind8                	popA_ind7                	0.081349577339648454     
DIST(1)          	28           	8                  	0                  	popA_ind9                	popA_ind1                	0.090441910531874198     
DIST(1)          	29           	8                  	1                  	popA_ind9                	popA_ind2                	0.095867900816094934     
DIST(1)          	30           	8                  	2                  	popA_ind9                	popA_ind3                	0.058549180362201331     
DIST(1)          	31           	8                  	3                  	popA_ind9                	popA_ind4                	0.1026648719420763       
DIST(1)          	32           	8                  	4                  	popA_ind9                	popA_ind5                	0.069954678492391426     
DIST(1)          	33           	8                  	5                  	popA_ind9                	popA_ind6                	0.092754438818220694     
DIST(1)          	34           	8                  	6                  	popA_ind9                	popA_ind7                	0.067255283583122355     
DIST(1)          	35           	8                  	7                  	popA_ind9                	popA_ind8                	0.081902410789416333     
DIST(1)          	36           	9                  	0                  	popA_ind10               	popA_ind1                	0.12644992674147423      
DIST(1)          	37           	9                  	1                  	popA_ind10               	popA_ind2                	0.0966360124762633       
DIST(1)          	38           	9                  	2                  	popA_ind10               	popA_ind3                	0.10716580364420696      
DIST(1)          	39           	9                  	3                  	popA_ind10               	popA_ind4                	0.101843573924255        
DIST(1)          	40           	9                  	4                  	popA_ind10               	popA_ind5                	0.10530742687306807      
DIST(1)          	41           	9                  	5                  	popA_ind10               	popA_ind6                	0.079663031092544578     
DIST(1)          	42           	9                  	6                  	popA_ind10               	popA_ind7                	0.13024246593697633      
DIST(1)          	43           	9                  	7                  	popA_ind10               	popA_ind8                	0.086512330866314727     
DIST(1)          	44           	9                  	8                  	popA_ind10               	popA_ind9                	0.091715136836088729     
DIST(1)          	45           	10                 	0                  	popB_ind1                	popA_ind1                	0.2336180839991483       
DIST(1)          	46           	10                 	1                  	popB_ind1                	popA_ind2                	0.18954759077640648      
DIST(1)          	47           	10                 	2                  	popB_ind1                	popA_ind3                	0.20283944424946287      
DIST(1)          	48           	10                 	3                  	popB_ind1                	popA_ind4                	0.23700842812959011      
DIST(1)          	49           	10                 	4                  	popB_ind1                	popA_ind5                	0.21720537528782002      
DIST(1)          	50           	10                 	5                  	popB_ind1                	popA_ind6                	0.16646532629003263      
DIST(1)          	51           	10                 	6                  	popB_ind1                	popA_ind7                	0.17940566801556734      
DIST(1)          	52           	10                 	7                  	popB_ind1                	popA_ind8                	0.18782183744554465      
DIST(1)          	53           	10                 	8                  	popB_ind1                	popA_ind9                	0.17880978453375065      
DIST(1)          	54           	10                 	9                  	popB_ind1                	popA_ind10               	0.1331576985459951       
DIST(1)          	55           	11                 	0                  	popB_ind2                	popA_ind1                	0.1557493052135199       
DIST(1)          	56           	11                 	1                  	popB_ind2                	popA_ind2                	0.12902711380341356      
DIST(1)          	57           	11                 	2                  	popB_ind2                	popA_ind3                	0.18546996912365249      
DIST(1)          	58           	11                 	3                  	popB_ind2                	popA_ind4                	0.20095151449527437      
DIST(1)          	59           	11                 	4                  	popB_ind2                	popA_ind5                	0.13921199302111537      
DIST(1)          	60           	11                 	5                  	popB_ind2                	popA_ind6                	0.11863830608516811      
DIST(1)          	61           	11                 	6                  	popB_ind2                	popA_ind7                	0.21865482444832432      
DIST(1)          	62           	11                 	7                  	popB_ind2                	popA_ind8                	0.16953188761612992      
DIST(1)          	63           	11                 	8                  	popB_ind2                	popA_ind9                	0.13355896813518423      
DIST(1)          	64           	11                 	9                  	popB_ind2                	popA_ind10               	0.13685061378832802      
DIST(1)          	65           	11                 	10                 	popB_ind2                	popB_ind1                	0.082160885097566277     
DIST(1)          	66           	12                 	0                  	popB_ind3                	popA_ind1                	0.18435734223985234      
DIST(1)          	67           	12                 	1                  	popB_ind3                	popA_ind2                	0.18073172098884416      
DIST(1)          	68           	12                 	2                  	popB_ind3                	popA_ind3                	0.12921894115135346      
DIST(1)          	69           	12                 	3                  	popB_ind3                	popA_ind4                	0.1228474585695065       
DIST(1)          	70           	12                 	4                  	popB_ind3                	popA_ind5                	0.18104606093364703      
DIST(1)          	71           	12                 	5                  	popB_ind3                	popA_ind6                	0.19471725344017074      
DIST(1)          	72           	12                 	6                  	popB_ind3                	popA_ind7                	0.1936455547548524       
DIST(1)          	73           	12                 	7                  	popB_ind3                	popA_ind8                	0.11517994569806499      
DIST(1)          	74           	12                 	8                  	popB_ind3                	popA_ind9                	0.1348604073845501       
DIST(1)          	75           	12                 	9                  	popB_ind3                	popA_ind10               	0.19398220480765643      
DIST(1)          	76           	12                 	10                 	popB_ind3                	popB_ind1                	0.16274535265075502      
DIST(1)          	77           	12                 	11                 	popB_ind3                	popB_ind2                	0.13373797837552709      
DIST(1)          	78           	13                 	0                  	popB_ind4                	popA_ind1                	0.13365602036443922      
DIST(1)          	79           	13                 	1                  	popB_ind4                	popA_ind2                	0.17171434159611618      
DIST(1)          	80           	13                 	2                  	popB_ind4                	popA_ind3                	0.14385416636385351      
DIST(1)          	81           	13                 	3                  	popB_ind4                	popA_ind4                	0.1362039295122727       
DIST(1)          	82           	13                 	4                  	popB_ind4                	popA_ind5                	0.13844854551032734      
DIST(1)          	83           	13                 	5                  	popB_ind4                	popA_ind6                	0.20507900555218048      
DIST(1)          	84           	13                 	6                  	popB_ind4                	popA_ind7                	0.16674800759215763      
DIST(1)          	85           	13                 	7                  	popB_ind4                	popA_ind8                	0.19158365508898345      
DIST(1)          	86           	13                 	8                  	popB_ind4                	popA_ind9                	0.14481723658856768      
DIST(1)          	87           	13                 	9                  	popB_ind4                	popA_ind10               	0.17216895148979977      
DIST(1)          	88           	13                 	10                 	popB_ind4                	popB_ind1                	0.10146661607076184      
DIST(1)          	89           	13                 	11                 	popB_ind4                	popB_ind2                	0.12618862661658065      
DIST(1)          	90           	13                 	12                 	popB_ind4                	popB_ind3                	0.058739516198463106     
DIST(1)          	91           	14                 	0                  	popB_ind5                	popA_ind1                	0.21820072728406747      
DIST(1)          	92           	14                 	1                  	popB_ind5                	popA_ind2                	0.14650675561445772      
DIST(1)          	93           	14                 	2                  	popB_ind5                	popA_ind3                	0.19196571222817471      
DIST(1)          	94           	14                 	3                  	popB_ind5                	popA_ind4                	0.17618931611687783      
DIST(1)          	95           	14                 	4                  	popB_ind5                	popA_ind5                	0.16128517489385924      
DIST(1)          	96           	14                 	5                  	popB_ind5                	popA_ind6                	0.15657767928982225      
DIST(1)          	97           	14                 	6                  	popB_ind5                	popA_ind7                	0.17567679918954107      
DIST(1)          	98           	14                 	7                  	popB_ind5                	popA_ind8                	0.16873546816436644      
DIST(1)          	99           	14                 	8                  	popB_ind5                	popA_ind9                	0.19334619577350504      
DIST(1)          	100          	14                 	9                  	popB_ind5                	popA_ind10               	0.14622659497598672      
DIST(1)          	101          	14                 	10                 	popB_ind5                	popB_ind1                	0.076306207154471964     
DIST(1)          	102          	14                 	11                 	popB_ind5                	popB_ind2                	0.066529877811036503     
DIST(1)          	103          	14                 	12                 	popB_ind5                	popB_ind3                	0.14094098427547055      
DIST(1)          	104          	14                 	13                 	popB_ind5                	popB_ind4                	0.096059387487813036     
DIST(1)          	105          	15                 	0                  	popB_ind6                	popA_ind1                	0.1384893437417391       
DIST(1)          	106          	15                 	1                  	popB_ind6                	popA_ind2                	0.15612309594158832      
DIST(1)          	107          	15                 	2                  	popB_ind6                	popA_ind3                	0.17415526086063793      
DIST(1)          	108          	15                 	3                  	popB_ind6                	popA_ind4                	0.19717244117112803      
DIST(1)          	109          	15                 	4                  	popB_ind6                	popA_ind5                	0.2018992075184749       
DIST(1)          	110          	15                 	5                  	popB_ind6                	popA_ind6                	0.21586937475895562      
DIST(1)          	111          	15                 	6                  	popB_ind6                	popA_ind7                	0.18762449822932747      
DIST(1)          	112          	15                 	7                  	popB_ind6                	popA_ind8                	0.19803684293859389      
DIST(1)          	113          	15                 	8                  	popB_ind6                	popA_ind9                	0.17157542180045365      
DIST(1)          	114          	15                 	9                  	popB_ind6                	popA_ind10               	0.15205420861488694      
DIST(1)          	115          	15                 	10                 	popB_ind6                	popB_ind1                	0.13700903211250814      
DIST(1)          	116          	15                 	11                 	popB_ind6                	popB_ind2                	0.08776638529312078      
DIST(1)          	117          	15                 	12                 	popB_ind6                	popB_ind3                	0.14859284705566389      
DIST(1)          	118          	15                 	13                 	popB_ind6                	popB_ind4                	0.13156715752246448      
DIST(1)          	119          	15                 	14                 	popB_ind6                	popB_ind5                	0.10742285964688825      
DIST(1)          	120          	16                 	0                  	popB_ind7                	popA_ind1                	0.17973240429294551      
DIST(1)          	121          	16                 	1                  	popB_ind7                	popA_ind2                	0.14078806644125305      
DIST(1)          	122          	16                 	2                  	popB_ind7                	popA_ind3                	0.091761684712605388     
DIST(1)          	123          	16                 	3                  	popB_ind7                	popA_ind4                	0.14915398801606841      
DIST(1)          	124          	16                 	4                  	popB_ind7                	popA_ind5                	0.16024652644638937      
DIST(1)          	125          	16                 	5                  	popB_ind7                	popA_ind6                	0.16262642621955042      
DIST(1)          	126          	16                 	6                  	popB_ind7                	popA_ind7                	0.1238278896798205       
DIST(1)          	127          	16                 	7                  	popB_ind7                	popA_ind8                	0.12664600605320162      
DIST(1)          	128          	16                 	8                  	popB_ind7                	popA_ind9                	0.11738137106945715      
DIST(1)          	129          	16                 	9                  	popB_ind7                	popA_ind10               	0.15210423771205944      
DIST(1)          	130          	16                 	10                 	popB_ind7                	popB_ind1                	0.082060786862680424     
DIST(1)          	131          	16                 	11                 	popB_ind7                	popB_ind2                	0.13608893161450789      
DIST(1)          	132          	16                 	12                 	popB_ind7                	popB_ind3                	0.080284930300773313     
DIST(1)          	133          	16                 	13                 	popB_ind7                	popB_ind4                	0.069619761483182191     
DIST(1)          	134          	16                 	14                 	popB_ind7                	popB_ind5                	0.13062814343793538      
DIST(1)          	135          	16                 	15                 	popB_ind7                	popB_ind6                	0.14559357684436566      
DIST(1)          	136          	17                 	0                  	popB_ind8                	popA_ind1                	0.16567511205365207      
DIST(1)          	137          	17                 	1                  	popB_ind8                	popA_ind2                	0.12640773566882468      
DIST(1)          	138          	17                 	2                  	popB_ind8                	popA_ind3                	0.17088661897786564      
DIST(1)          	139          	17                 	3                  	popB_ind8                	popA_ind4                	0.23133671988418381      
DIST(1)          	140          	17                 	4                  	popB_ind8                	popA_ind5                	0.19225338752473931      
DIST(1)          	141          	17                 	5                  	popB_ind8                	popA_ind6                	0.19726913409919913      
DIST(1)          	142          	17                 	6                  	popB_ind8                	popA_ind7                	0.14794873930210076      
DIST(1)          	143          	17                 	7                  	popB_ind8                	popA_ind8                	0.17227982719132295      
DIST(1)          	144          	17                 	8                  	popB_ind8                	popA_ind9                	0.14959957714123034      
DIST(1)          	145          	17                 	9                  	popB_ind8                	popA_ind10               	0.15249983598143488      
DIST(1)          	146          	17                 	10                 	popB_ind8                	popB_ind1                	0.096185128202881956     
DIST(1)          	147          	17                 	11                 	popB_ind8                	popB_ind2                	0.12356963919485967      
DIST(1)          	148          	17                 	12                 	popB_ind8                	popB_ind3                	0.097396995678103074     
DIST(1)          	149          	17                 	13                 	popB_ind8                	popB_ind4                	0.11707025754164227      
DIST(1)          	150          	17                 	14                 	popB_ind8                	popB_ind5                	0.077202720953757867     
DIST(1)          	151          	17                 	15                 	popB_ind8                	popB_ind6                	0.12856793997731336      
DIST(1)          	152          	17                 	16                 	popB_ind8                	popB_ind7                	0.12011526516704352      
DIST(1)          	153          	18                 	0                  	popB_ind9                	popA_ind1                	0.12538366546380877      
DIST(1)          	154          	18                 	1                  	popB_ind9                	popA_ind2                	0.15367870315230986      
DIST(1)          	155          	18                 	2                  	popB_ind9                	popA_ind3                	0.10515210094250806      
DIST(1)          	156          	18                 	3                  	popB_ind9                	popA_ind4                	0.15781540012958445      
DIST(1)          	157          	18                 	4                  	popB_ind9                	popA_ind5                	0.12211001480676478      
DIST(1)          	158          	18                 	5                  	popB_ind9                	popA_ind6                	0.13742681971106563      
DIST(1)          	159          	18                 	6                  	popB_ind9                	popA_ind7                	0.11147969497017537      
DIST(1)          	160          	18                 	7                  	popB_ind9                	popA_ind8                	0.14821969503273408      
DIST(1)          	161          	18                 	8                  	popB_ind9                	popA_ind9                	0.095865949319817359     
DIST(1)          	162          	18                 	9                  	popB_ind9                	popA_ind10               	0.15796362529584293      
DIST(1)          	163          	18                 	10                 	popB_ind9                	popB_ind1                	0.11159237196376755      
DIST(1)          	164          	18                 	11                 	popB_ind9                	popB_ind2                	0.10566166338278882      
DIST(1)          	165          	18                 	12                 	popB_ind9                	popB_ind3                	0.11834799348704117      
DIST(1)          	166          	18                 	13                 	popB_ind9                	popB_ind4                	0.09777857849635356      
DIST(1)          	167          	18                 	14                 	popB_ind9                	popB_ind5                	0.092115999449317851     
DIST(1)          	168          	18                 	15                 	popB_ind9                	popB_ind6                	0.11123602906816323      
DIST(1)          	169          	18                 	16                 	popB_ind9                	popB_ind7                	0.079927548958855135     
DIST(1)          	170          	18                 	17                 	popB_ind9                	popB_ind8                	0.097300091716630593     
DIST(1)          	171          	19                 	0                  	popB_ind10               	popA_ind1                	0.21307666193115798      
DIST(1)          	172          	19                 	1                  	popB_ind10               	popA_ind2                	0.1703770676526504       
DIST(1)          	173          	19                 	2                  	popB_ind10               	popA_ind3                	0.14639700483239174      
DIST(1)          	174          	19                 	3                  	popB_ind10               	popA_ind4                	0.19214745347349027      
DIST(1)          	175          	19                 	4                  	popB_ind10               	popA_ind5                	0.13061699069084554      
DIST(1)          	176          	19                 	5                  	popB_ind10               	popA_ind6                	0.15657229526571761      
DIST(1)          	177          	19                 	6                  	popB_ind10               	popA_ind7                	0.16985562811883512      
DIST(1)          	178          	19                 	7                  	popB_ind10               	popA_ind8                	0.18171893336639938      
DIST(1)          	179          	19                 	8                  	popB_ind10               	popA_ind9                	0.1671966477195756       
DIST(1)          	180          	19                 	9                  	popB_ind10               	popA_ind10               	0.16347371290641377      
DIST(1)          	181          	19                 	10                 	popB_ind10               	popB_ind1                	0.12687343944407431      
DIST(1)          	182          	19                 	11                 	popB_ind10               	popB_ind2                	0.14430340760295579      
DIST(1)          	183          	19                 	12                 	popB_ind10               	popB_ind3                	0.071379753531828252     
DIST(1)          	184          	19                 	13                 	popB_ind10               	popB_ind4                	0.074913241182063356     
DIST(1)          	185          	19                 	14                 	popB_ind10               	popB_ind5                	0.11966181426920416      
DIST(1)          	186          	19                 	15                 	popB_ind10               	popB_ind6                	0.14635490020317699      
DIST(1)          	187          	19                 	16                 	popB_ind10               	popB_ind7                	0.093136998894170472     
DIST(1)          	188          	19                 	17                 	popB_ind10               	popB_ind8                	0.079366584664830125     
DIST(1)          	189          	19                 	18                 	popB_ind10               	popB_ind9                	0.11427538896625043      
DIST(1)          	190          	20                 	0                  	popC_ind1                	popA_ind1                	0.2015332824908172       
DIST(1)          	191          	20                 	1                  	popC_ind1                	popA_ind2                	0.21319717169203276      
DIST(1)          	192          	20                 	2                  	popC_ind1                	popA_ind3                	0.19935098429448844      
DIST(1)          	193          	20                 	3                  	popC_ind1                	popA_ind4                	0.22274743747405051      
DIST(1)          	194          	20                 	4                  	popC_ind1                	popA_ind5                	0.13482631394144035      
DIST(1)          	195          	20                 	5                  	popC_ind1                	popA_ind6                	0.22349667533920595      
DIST(1)          	196          	20                 	6                  	popC_ind1                	popA_ind7                	0.21226372044518796      
DIST(1)          	197          	20                 	7                  	popC_ind1                	popA_ind8                	0.15416149791707118      
DIST(1)          	198          	20                 	8                  	popC_ind1                	popA_ind9                	0.14624866514036183      
DIST(1)          	199          	20                 	9                  	popC_ind1                	popA_ind10               	0.19273603119366223      
DIST(1)          	200          	20                 	10                 	popC_ind1                	popB_ind1                	0.18992584418755426      
DIST(1)          	201          	20                 	11                 	popC_ind1                	popB_ind2                	0.13916576111370926      
DIST(1)          	202          	20                 	12                 	popC_ind1                	popB_ind3                	0.21871802943394594      
DIST(1)          	203          	20                 	13                 	popC_ind1                	popB_ind4                	0.23894372190687277      
DIST(1)          	204          	20                 	14                 	popC_ind1                	popB_ind5                	0.17615682862512969      
DIST(1)          	205          	20                 	15                 	popC_ind1                	popB_ind6                	0.21229738138034671      
DIST(1)          	206          	20                 	16                 	popC_ind1                	popB_ind7                	0.25276437247717154      
DIST(1)          	207          	20                 	17                 	popC_ind1                	popB_ind8                	0.22832392155407683      
DIST(1)          	208          	20                 	18                 	popC_ind1                	popB_ind9                	0.20775044417262878      
DIST(1)          	209          	20                 	19                 	popC_ind1                	popB_ind10               	0.14419441349830023      
DIST(1)          	210          	21                 	0                  	popC_ind2                	popA_ind1                	0.22010928878735497      
DIST(1)          	211          	21                 	1                  	popC_ind2                	popA_ind2                	0.19596048314884862      
DIST(1)          	212          	21                 	2                  	popC_ind2                	popA_ind3                	0.13348044864677666      
DIST(1)          	213          	21                 	3                  	popC_ind2                	popA_ind4                	0.22768146183832752      
DIST(1)          	214          	21                 	4                  	popC_ind2                	popA_ind5                	0.17210693161572238      
DIST(1)          	215          	21                 	5                  	popC_ind2                	popA_ind6                	0.22974283767662584      
DIST(1)          	216          	21                 	6                  	popC_ind2                	popA_ind7                	0.20868695076057764      
DIST(1)          	217          	21                 	7                  	popC_ind2                	popA_ind8                	0.21980391670649738      
DIST(1)          	218          	21                 	8                  	popC_ind2                	popA_ind9                	0.20436847500812061      
DIST(1)          	219          	21                 	9                  	popC_ind2                	popA_ind10               	0.2377392532730476       
DIST(1)          	220          	21                 	10                 	popC_ind2                	popB_ind1                	0.24110193894173057      
DIST(1)          	221          	21                 	11                 	popC_ind2                	popB_ind2                	0.18274924278785168      
DIST(1)          	222          	21                 	12                 	popC_ind2                	popB_ind3                	0.22482126011468134      
DIST(1)          	223          	21                 	13                 	popC_ind2                	popB_ind4                	0.21303988719470013      
DIST(1)          	224          	21                 	14                 	popC_ind2                	popB_ind5                	0.24704618511131377      
DIST(1)          	225          	21                 	15                 	popC_ind2                	popB_ind6                	0.17241113608750086      
DIST(1)          	226          	21                 	16                 	popC_ind2                	popB_ind7                	0.22743761668678641      
DIST(1)          	227          	21                 	17                 	popC_ind2                	popB_ind8                	0.24500901223963265      
DIST(1)          	228          	21                 	18                 	popC_ind2                	popB_ind9                	0.17401907242445566      
DIST(1)          	229          	21                 	19                 	popC_ind2                	popB_ind10               	0.17641728721519209      
DIST(1)          	230          	21                 	20                 	popC_ind2                	popC_ind1                	0.085929952914894581     
DIST(1)          	231          	22                 	0                  	popC_ind3                	popA_ind1                	0.20992072881488483      
DIST(1)          	232          	22                 	1                  	popC_ind3                	popA_ind2                	0.19856025357945117      
DIST(1)          	233          	22                 	2                  	popC_ind3                	popA_ind3                	0.18137186054747734      
DIST(1)          	234          	22                 	3                  	popC_ind3                	popA_ind4                	0.22563648648046841      
DIST(1)          	235          	22                 	4                  	popC_ind3                	popA_ind5                	0.19355061640016896      
DIST(1)          	236          	22                 	5                  	popC_ind3                	popA_ind6                	0.20568966173104594      
DIST(1)          	237          	22                 	6                  	popC_ind3                	popA_ind7                	0.18097502863792991      
DIST(1)          	238          	22                 	7                  	popC_ind3                	popA_ind8                	0.19008316661633073      
DIST(1)          	239          	22                 	8                  	popC_ind3                	popA_ind9                	0.20914506792198156      
DIST(1)          	240          	22                 	9                  	popC_ind3                	popA_ind10               	0.19661696569243406      
DIST(1)          	241          	22                 	10                 	popC_ind3                	popB_ind1                	0.20516226119894762      
DIST(1)          	242          	22                 	11                 	popC_ind3                	popB_ind2                	0.17419158370645485      
DIST(1)          	243          	22                 	12                 	popC_ind3                	popB_ind3                	0.223667855964443        
DIST(1)          	244          	22                 	13                 	popC_ind3                	popB_ind4                	0.23600090782376576      
DIST(1)          	245          	22                 	14                 	popC_ind3                	popB_ind5                	0.24198609467351978      
DIST(1)          	246          	22                 	15                 	popC_ind3                	popB_ind6                	0.19888742525337749      
DIST(1)          	247          	22                 	16                 	popC_ind3                	popB_ind7                	0.21567635282876937      
DIST(1)          	248          	22                 	17                 	popC_ind3                	popB_ind8                	0.20019062125200154      
DIST(1)          	249          	22                 	18                 	popC_ind3                	popB_ind9                	0.20930096488445007      
DIST(1)          	250          	22                 	19                 	popC_ind3                	popB_ind10               	0.19961998820517379      
DIST(1)          	251          	22                 	20                 	popC_ind3                	popC_ind1                	0.06458724038052234      
DIST(1)          	252          	22                 	21                 	popC_ind3                	popC_ind2                	0.083556136378852577     
DIST(1)          	253          	23                 	0                  	popC_ind4                	popA_ind1                	0.18011973954276245      
DIST(1)          	254          	23                 	1                  	popC_ind4                	popA_ind2                	0.17799342432036946      
DIST(1)          	255          	23                 	2                  	popC_ind4                	popA_ind3                	0.14461749966416013      
DIST(1)          	256          	23                 	3                  	popC_ind4                	popA_ind4                	0.19707754359541055      
DIST(1)          	257          	23                 	4                  	popC_ind4                	popA_ind5                	0.081157311261563378     
DIST(1)          	258          	23                 	5                  	popC_ind4                	popA_ind6                	0.19920225108790207      
DIST(1)          	259          	23                 	6                  	popC_ind4                	popA_ind7                	0.23824256705359412      
DIST(1)          	260          	23                 	7                  	popC_ind4                	popA_ind8                	0.12737489517049941      
DIST(1)          	261          	23                 	8                  	popC_ind4                	popA_ind9                	0.13183984482082806      
DIST(1)          	262          	23                 	9                  	popC_ind4                	popA_ind10               	0.16902390438926712      
DIST(1)          	263          	23                 	10                 	popC_ind4                	popB_ind1                	0.16848696867103963      
DIST(1)          	264          	23                 	11                 	popC_ind4                	popB_ind2                	0.1485571300127603       
DIST(1)          	265          	23                 	12                 	popC_ind4                	popB_ind3                	0.22618623471090815      
DIST(1)          	266          	23                 	13                 	popC_ind4                	popB_ind4                	0.19733227421567021      
DIST(1)          	267          	23                 	14                 	popC_ind4                	popB_ind5                	0.21928494446825603      
DIST(1)          	268          	23                 	15                 	popC_ind4                	popB_ind6                	0.14381946282788433      
DIST(1)          	269          	23                 	16                 	popC_ind4                	popB_ind7                	0.18400919938333196      
DIST(1)          	270          	23                 	17                 	popC_ind4                	popB_ind8                	0.2153539390760478       
DIST(1)          	271          	23                 	18                 	popC_ind4                	popB_ind9                	0.17428643282842946      
DIST(1)          	272          	23                 	19                 	popC_ind4                	popB_ind10               	0.1677931257830779       
DIST(1)          	273          	23                 	20                 	popC_ind4                	popC_ind1                	0.088996394857584066     
DIST(1)          	274          	23                 	21                 	popC_ind4                	popC_ind2                	0.0636691408407625       
DIST(1)          	275          	23                 	22                 	popC_ind4                	popC_ind3                	0.087365595598285314     
DIST(1)          	276          	24                 	0                  	popC_ind5                	popA_ind1                	0.20325323130980855      
DIST(1)          	277          	24                 	1                  	popC_ind5                	popA_ind2                	0.14702450642836418      
DIST(1)          	278          	24                 	2                  	popC_ind5                	popA_ind3                	0.12378919856248129      
DIST(1)          	279          	24                 	3                  	popC_ind5                	popA_ind4                	0.16490775905017879      
DIST(1)          	280          	24                 	4                  	popC_ind5                	popA_ind5                	0.12756045244004607      
DIST(1)          	281          	24                 	5                  	popC_ind5                	popA_ind6                	0.15275739577477193      
DIST(1)          	282          	24                 	6                  	popC_ind5                	popA_ind7                	0.14767875278433951      
DIST(1)          	283          	24                 	7                  	popC_ind5                	popA_ind8                	0.11134245564017298      
DIST(1)          	284          	24                 	8                  	popC_ind5                	popA_ind9                	0.1545797465156945       
DIST(1)          	285          	24                 	9                  	popC_ind5                	popA_ind10               	0.16103514198867441      
DIST(1)          	286          	24                 	10                 	popC_ind5                	popB_ind1                	0.18153955539439667      
DIST(1)          	287          	24                 	11                 	popC_ind5                	popB_ind2                	0.12645621735530049      
DIST(1)          	288          	24                 	12                 	popC_ind5                	popB_ind3                	0.21219604150214644      
DIST(1)          	289          	24                 	13                 	popC_ind5                	popB_ind4                	0.22434630118986323      
DIST(1)          	290          	24                 	14                 	popC_ind5                	popB_ind5                	0.21254986108129204      
DIST(1)          	291          	24                 	15                 	popC_ind5                	popB_ind6                	0.16881940483549152      
DIST(1)          	292          	24                 	16                 	popC_ind5                	popB_ind7                	0.14723056087689204      
DIST(1)          	293          	24                 	17                 	popC_ind5                	popB_ind8                	0.20735358341427781      
DIST(1)          	294          	24                 	18                 	popC_ind5                	popB_ind9                	0.15812908243462728      
DIST(1)          	295          	24                 	19                 	popC_ind5                	popB_ind10               	0.12767489028796097      
DIST(1)          	296          	24                 	20                 	popC_ind5                	popC_ind1                	0.038373039060623483     
DIST(1)          	297          	24                 	21                 	popC_ind5                	popC_ind2                	0.092059644589807504     
DIST(1)          	298          	24                 	22                 	popC_ind5                	popC_ind3                	0.074413818107818774     
DIST(1)          	299          	24                 	23                 	popC_ind5                	popC_ind4                	0.053726777929373526     
DIST(1)          	300          	25                 	0                  	popC_ind6                	popA_ind1                	0.20984031673621628      
DIST(1)          	301          	25                 	1                  	popC_ind6                	popA_ind2                	0.17442468605602213      
DIST(1)          	302          	25                 	2                  	popC_ind6                	popA_ind3                	0.25895711393091969      
DIST(1)          	303          	25                 	3                  	popC_ind6                	popA_ind4                	0.2116313684226776       
DIST(1)          	304          	25                 	4                  	popC_ind6                	popA_ind5                	0.16355731741209173      
DIST(1)          	305          	25                 	5                  	popC_ind6                	popA_ind6                	0.26326882256131307      
DIST(1)          	306          	25                 	6                  	popC_ind6                	popA_ind7                	0.26393090198139951      
DIST(1)          	307          	25                 	7                  	popC_ind6                	popA_ind8                	0.168157258263574        
DIST(1)          	308          	25                 	8                  	popC_ind6                	popA_ind9                	0.22714663644873589      
DIST(1)          	309          	25                 	9                  	popC_ind6                	popA_ind10               	0.24610165008253934      
DIST(1)          	310          	25                 	10                 	popC_ind6                	popB_ind1                	0.24285448097296566      
DIST(1)          	311          	25                 	11                 	popC_ind6                	popB_ind2                	0.18794323710484434      
DIST(1)          	312          	25                 	12                 	popC_ind6                	popB_ind3                	0.2442810931287237       
DIST(1)          	313          	25                 	13                 	popC_ind6                	popB_ind4                	0.18087135708310473      
DIST(1)          	314          	25                 	14                 	popC_ind6                	popB_ind5                	0.18653439343397099      
DIST(1)          	315          	25                 	15                 	popC_ind6                	popB_ind6                	0.2058210967759192       
DIST(1)          	316          	25                 	16                 	popC_ind6                	popB_ind7                	0.22320203943759701      
DIST(1)          	317          	25                 	17                 	popC_ind6                	popB_ind8                	0.21083984955048851      
DIST(1)          	318          	25                 	18                 	popC_ind6                	popB_ind9                	0.18504443247286231      
DIST(1)          	319          	25                 	19                 	popC_ind6                	popB_ind10               	0.1909788577390174       
DIST(1)          	320          	25                 	20                 	popC_ind6                	popC_ind1                	0.061336527554696978     
DIST(1)          	321          	25                 	21                 	popC_ind6                	popC_ind2                	0.060911166662706975     
DIST(1)          	322          	25                 	22                 	popC_ind6                	popC_ind3                	0.088665456548880514     
DIST(1)          	323          	25                 	23                 	popC_ind6                	popC_ind4                	0.07408113418657028      
DIST(1)          	324          	25                 	24                 	popC_ind6                	popC_ind5                	0.053431671951617882     
DIST(1)          	325          	26                 	0                  	popC_ind7                	popA_ind1                	0.23655857744765849      
DIST(1)          	326          	26                 	1                  	popC_ind7                	popA_ind2                	0.24830740085268679      
DIST(1)          	327          	26                 	2                  	popC_ind7                	popA_ind3                	0.22165859931328433      
DIST(1)          	328          	26                 	3                  	popC_ind7                	popA_ind4                	0.18967820697965587      
DIST(1)          	329          	26                 	4                  	popC_ind7                	popA_ind5                	0.20138998552086085      
DIST(1)          	330          	26                 	5                  	popC_ind7                	popA_ind6                	0.28795560981523416      
DIST(1)          	331          	26                 	6                  	popC_ind7                	popA_ind7                	0.26873541148604974      
DIST(1)          	332          	26                 	7                  	popC_ind7                	popA_ind8                	0.19675588091388402      
DIST(1)          	333          	26                 	8                  	popC_ind7                	popA_ind9                	0.22811461996768165      
DIST(1)          	334          	26                 	9                  	popC_ind7                	popA_ind10               	0.20867271024605341      
DIST(1)          	335          	26                 	10                 	popC_ind7                	popB_ind1                	0.24061109001033376      
DIST(1)          	336          	26                 	11                 	popC_ind7                	popB_ind2                	0.19877923217784804      
DIST(1)          	337          	26                 	12                 	popC_ind7                	popB_ind3                	0.21547008840665588      
DIST(1)          	338          	26                 	13                 	popC_ind7                	popB_ind4                	0.19767769172169303      
DIST(1)          	339          	26                 	14                 	popC_ind7                	popB_ind5                	0.27509638561703653      
DIST(1)          	340          	26                 	15                 	popC_ind7                	popB_ind6                	0.25374463338126252      
DIST(1)          	341          	26                 	16                 	popC_ind7                	popB_ind7                	0.17857892695742442      
DIST(1)          	342          	26                 	17                 	popC_ind7                	popB_ind8                	0.22978234231581421      
DIST(1)          	343          	26                 	18                 	popC_ind7                	popB_ind9                	0.22425103724399476      
DIST(1)          	344          	26                 	19                 	popC_ind7                	popB_ind10               	0.19836430073578476      
DIST(1)          	345          	26                 	20                 	popC_ind7                	popC_ind1                	0.080468945192271074     
DIST(1)          	346          	26                 	21                 	popC_ind7                	popC_ind2                	0.077573278036358162     
DIST(1)          	347          	26                 	22                 	popC_ind7                	popC_ind3                	0.096624243274473581     
DIST(1)          	348          	26                 	23                 	popC_ind7                	popC_ind4                	0.078420141511245087     
DIST(1)          	349          	26                 	24                 	popC_ind7                	popC_ind5                	0.085692851086862698     
DIST(1)          	350          	26                 	25                 	popC_ind7                	popC_ind6                	0.094015717217159261     
DIST(1)          	351          	27                 	0                  	popC_ind8                	popA_ind1                	0.231717259781694        
DIST(1)          	352          	27                 	1                  	popC_ind8                	popA_ind2                	0.22090372507239445      
DIST(1)          	353          	27                 	2                  	popC_ind8                	popA_ind3                	0.2080977742654021       
DIST(1)          	354          	27                 	3                  	popC_ind8                	popA_ind4                	0.19897184497235881      
DIST(1)          	355          	27                 	4                  	popC_ind8                	popA_ind5                	0.20981466147567945      
DIST(1)          	356          	27                 	5                  	popC_ind8                	popA_ind6                	0.23941607366271042      
DIST(1)          	357          	27                 	6                  	popC_ind8                	popA_ind7                	0.20018477116636779      
DIST(1)          	358          	27                 	7                  	popC_ind8                	popA_ind8                	0.17347278125565413      
DIST(1)          	359          	27                 	8                  	popC_ind8                	popA_ind9                	0.20419846762542376      
DIST(1)          	360          	27                 	9                  	popC_ind8                	popA_ind10               	0.21781060486206255      
DIST(1)          	361          	27                 	10                 	popC_ind8                	popB_ind1                	0.22480874276002183      
DIST(1)          	362          	27                 	11                 	popC_ind8                	popB_ind2                	0.27709420638695609      
DIST(1)          	363          	27                 	12                 	popC_ind8                	popB_ind3                	0.22787476566118273      
DIST(1)          	364          	27                 	13                 	popC_ind8                	popB_ind4                	0.19083714046180264      
DIST(1)          	365          	27                 	14                 	popC_ind8                	popB_ind5                	0.30645922919112695      
DIST(1)          	366          	27                 	15                 	popC_ind8                	popB_ind6                	0.23059273142635173      
DIST(1)          	367          	27                 	16                 	popC_ind8                	popB_ind7                	0.23069574344760815      
DIST(1)          	368          	27                 	17                 	popC_ind8                	popB_ind8                	0.26584523918957531      
DIST(1)          	369          	27                 	18                 	popC_ind8                	popB_ind9                	0.22086612680901907      
DIST(1)          	370          	27                 	19                 	popC_ind8                	popB_ind10               	0.199701679793234        
DIST(1)          	371          	27                 	20                 	popC_ind8                	popC_ind1                	0.038367044363519689     
DIST(1)          	372          	27                 	21                 	popC_ind8                	popC_ind2                	0.057603914036350297     
DIST(1)          	373          	27                 	22                 	popC_ind8                	popC_ind3                	0.10002365391216403      
DIST(1)          	374          	27                 	23                 	popC_ind8                	popC_ind4                	0.074811407190029849     
DIST(1)          	375          	27                 	24                 	popC_ind8                	popC_ind5                	0.049541727827786659     
DIST(1)          	376          	27                 	25                 	popC_ind8                	popC_ind6                	0.083376674419933502     
DIST(1)          	377          	27                 	26                 	popC_ind8                	popC_ind7                	0.10231720437086897      
DIST(1)          	378          	28                 	0                  	popC_ind9                	popA_ind1                	0.27893889901517854      
DIST(1)          	379          	28                 	1                  	popC_ind9                	popA_ind2                	0.21410309795138516      
DIST(1)          	380          	28                 	2                  	popC_ind9                	popA_ind3                	0.19286296697013494      
DIST(1)          	381          	28                 	3                  	popC_ind9                	popA_ind4                	0.18686749308825443      
DIST(1)          	382          	28                 	4                  	popC_ind9                	popA_ind5                	0.15456903913146036      
DIST(1)          	383          	28                 	5                  	popC_ind9                	popA_ind6                	0.21061187153154848      
DIST(1)          	384          	28                 	6                  	popC_ind9                	popA_ind7                	0.22837868173009906      
DIST(1)          	385          	28                 	7                  	popC_ind9                	popA_ind8                	0.16378460291296684      
DIST(1)          	386          	28                 	8                  	popC_ind9                	popA_ind9                	0.17143383683045693      
DIST(1)          	387          	28                 	9                  	popC_ind9                	popA_ind10               	0.14333689465982724      
DIST(1)          	388          	28                 	10                 	popC_ind9                	popB_ind1                	0.24031571106900579      
DIST(1)          	389          	28                 	11                 	popC_ind9                	popB_ind2                	0.19753342149877429      
DIST(1)          	390          	28                 	12                 	popC_ind9                	popB_ind3                	0.24909517760238459      
DIST(1)          	391          	28                 	13                 	popC_ind9                	popB_ind4                	0.21178624032341026      
DIST(1)          	392          	28                 	14                 	popC_ind9                	popB_ind5                	0.23701875375463241      
DIST(1)          	393          	28                 	15                 	popC_ind9                	popB_ind6                	0.20842888691417302      
DIST(1)          	394          	28                 	16                 	popC_ind9                	popB_ind7                	0.27955586354207179      
DIST(1)          	395          	28                 	17                 	popC_ind9                	popB_ind8                	0.29242445900054814      
DIST(1)          	396          	28                 	18                 	popC_ind9                	popB_ind9                	0.25257668494528773      
DIST(1)          	397          	28                 	19                 	popC_ind9                	popB_ind10               	0.195070428112468        
DIST(1)          	398          	28                 	20                 	popC_ind9                	popC_ind1                	0.062541815271601547     
DIST(1)          	399          	28                 	21                 	popC_ind9                	popC_ind2                	0.096198048711564496     
DIST(1)          	400          	28                 	22                 	popC_ind9                	popC_ind3                	0.12801917742055624      
DIST(1)          	401          	28                 	23                 	popC_ind9                	popC_ind4                	0.095481720882076709     
DIST(1)          	402          	28                 	24                 	popC_ind9                	popC_ind5                	0.046700724870550984     
DIST(1)          	403          	28                 	25                 	popC_ind9                	popC_ind6                	0.078023295539278517     
DIST(1)          	404          	28                 	26                 	popC_ind9                	popC_ind7                	0.094692791475228286     
DIST(1)          	405          	28                 	27                 	popC_ind9                	popC_ind8                	0.094952695531240541     
DIST(1)          	406          	29                 	0                  	popC_ind10               	popA_ind1                	0.26359692965102527      
DIST(1)          	407          	29                 	1                  	popC_ind10               	popA_ind2                	0.24990291569905129      
DIST(1)          	408          	29                 	2                  	popC_ind10               	popA_ind3                	0.24506414677891714      
DIST(1)          	409          	29                 	3                  	popC_ind10               	popA_ind4                	0.27575413391200854      
DIST(1)          	410          	29                 	4                  	popC_ind10               	popA_ind5                	0.24426070118934101      
DIST(1)          	411          	29                 	5                  	popC_ind10               	popA_ind6                	0.22138250336297258      
DIST(1)          	412          	29                 	6                  	popC_ind10               	popA_ind7                	0.22809973495139563      
DIST(1)          	413          	29                 	7                  	popC_ind10               	popA_ind8                	0.22615354130733084      
DIST(1)          	414          	29                 	8                  	popC_ind10               	popA_ind9                	0.24653119815044344      
DIST(1)          	415          	29                 	9                  	popC_ind10               	popA_ind10               	0.22949346461278047      
DIST(1)          	416          	29                 	10                 	popC_ind10               	popB_ind1                	0.20338124312427924      
DIST(1)          	417          	29                 	11                 	popC_ind10               	popB_ind2                	0.23936921589388915      
DIST(1)          	418          	29                 	12                 	popC_ind10               	popB_ind3                	0.24991241331889913      
DIST(1)          	419          	29                 	13                 	popC_ind10               	popB_ind4                	0.22499187380586272      
DIST(1)          	420          	29                 	14                 	popC_ind10               	popB_ind5                	0.27313820987458215      
DIST(1)          	421          	29                 	15                 	popC_ind10               	popB_ind6                	0.25544194022353878      
DIST(1)          	422          	29                 	16                 	popC_ind10               	popB_ind7                	0.27546492200799055      
DIST(1)          	423          	29                 	17                 	popC_ind10               	popB_ind8                	0.26841693667885036      
DIST(1)          	424          	29                 	18                 	popC_ind10               	popB_ind9                	0.21351601714660889      
DIST(1)          	425          	29                 	19                 	popC_ind10               	popB_ind10               	0.19710585973442807      
DIST(1)          	426          	29                 	20                 	popC_ind10               	popC_ind1                	0.084133763987661134     
DIST(1)          	427          	29                 	21                 	popC_ind10               	popC_ind2                	0.089002603093127963     
DIST(1)          	428          	29                 	22                 	popC_ind10               	popC_ind3                	0.090197021640572708     
DIST(1)          	429          	29                 	23                 	popC_ind10               	popC_ind4                	0.070476984843123952     
DIST(1)          	430          	29                 	24                 	popC_ind10               	popC_ind5                	0.069405838498509634     
DIST(1)          	431          	29                 	25                 	popC_ind10               	popC_ind6                	0.081919128661087498     
DIST(1)          	432          	29                 	26                 	popC_ind10               	popC_ind7                	0.10398249075452973      
DIST(1)          	433          	29                 	27                 	popC_ind10               	popC_ind8                	0.088913855721502205     
DIST(1)          	434          	29                 	28                 	popC_ind10               	popC_ind9                	0.061781873583093866     
DIST(1)          	435          	30                 	0                  	popD_ind1                	popA_ind1                	0.15576874376073427      
DIST(1)          	436          	30                 	1                  	popD_ind1                	popA_ind2                	0.18354603047325474      
DIST(1)          	437          	30                 	2                  	popD_ind1                	popA_ind3                	0.15548244825433163      
DIST(1)          	438          	30                 	3                  	popD_ind1                	popA_ind4                	0.23636189166322452      
DIST(1)          	439          	30                 	4                  	popD_ind1                	popA_ind5                	0.16065407838892959      
DIST(1)          	440          	30                 	5                  	popD_ind1                	popA_ind6                	0.25146596833430707      
DIST(1)          	441          	30                 	6                  	popD_ind1                	popA_ind7                	0.18426987979109441      
DIST(1)          	442          	30                 	7                  	popD_ind1                	popA_ind8                	0.14242565452649741      
DIST(1)          	443          	30                 	8                  	popD_ind1                	popA_ind9                	0.18239425024363826      
DIST(1)          	444          	30                 	9                  	popD_ind1                	popA_ind10               	0.20453525675933965      
DIST(1)          	445          	30                 	10                 	popD_ind1                	popB_ind1                	0.20477225865304255      
DIST(1)          	446          	30                 	11                 	popD_ind1                	popB_ind2                	0.20888693041062803      
DIST(1)          	447          	30                 	12                 	popD_ind1                	popB_ind3                	0.19351193017586943      
DIST(1)          	448          	30                 	13                 	popD_ind1                	popB_ind4                	0.19991354311571102      
DIST(1)          	449          	30                 	14                 	popD_ind1                	popB_ind5                	0.24126918206266984      
DIST(1)          	450          	30                 	15                 	popD_ind1                	popB_ind6                	0.17394113001215425      
DIST(1)          	451          	30                 	16                 	popD_ind1                	popB_ind7                	0.18805765501821284      
DIST(1)          	452          	30                 	17                 	popD_ind1                	popB_ind8                	0.18541141583156168      
DIST(1)          	453          	30                 	18                 	popD_ind1                	popB_ind9                	0.15971641445464069      
DIST(1)          	454          	30                 	19                 	popD_ind1                	popB_ind10               	0.20158372573305069      
DIST(1)          	455          	30                 	20                 	popD_ind1                	popC_ind1                	0.12984999798110181      
DIST(1)          	456          	30                 	21                 	popD_ind1                	popC_ind2                	0.093924844356783097     
DIST(1)          	457          	30                 	22                 	popD_ind1                	popC_ind3                	0.1495852717103146       
DIST(1)          	458          	30                 	23                 	popD_ind1                	popC_ind4                	0.087817162380279229     
DIST(1)          	459          	30                 	24                 	popD_ind1                	popC_ind5                	0.074551563154910003     
DIST(1)          	460          	30                 	25                 	popD_ind1                	popC_ind6                	0.11993332403573874      
DIST(1)          	461          	30                 	26                 	popD_ind1                	popC_ind7                	0.16517261220996163      
DIST(1)          	462          	30                 	27                 	popD_ind1                	popC_ind8                	0.11194600643558583      
DIST(1)          	463          	30                 	28                 	popD_ind1                	popC_ind9                	0.13796409469097354      
DIST(1)          	464          	30                 	29                 	popD_ind1                	popC_ind10               	0.16068347267613489      
DIST(1)          	465          	31                 	0                  	popD_ind2                	popA_ind1                	0.2252919853265799       
DIST(1)          	466          	31                 	1                  	popD_ind2                	popA_ind2                	0.15273144787650705      
DIST(1)          	467          	31                 	2                  	popD_ind2                	popA_ind3                	0.23137037387074694      
DIST(1)          	468          	31                 	3                  	popD_ind2                	popA_ind4                	0.21336392171997709      
DIST(1)          	469          	31                 	4                  	popD_ind2                	popA_ind5                	0.18413506961360537      
DIST(1)          	470          	31                 	5                  	popD_ind2                	popA_ind6                	0.20093618363776039      
DIST(1)          	471          	31                 	6                  	popD_ind2                	popA_ind7                	0.17424317963019309      
DIST(1)          	472          	31                 	7                  	popD_ind2                	popA_ind8                	0.22066812525444188      
DIST(1)          	473          	31                 	8                  	popD_ind2                	popA_ind9                	0.16343172198789696      
DIST(1)          	474          	31                 	9                  	popD_ind2                	popA_ind10               	0.26837875257604882      
DIST(1)          	475          	31                 	10                 	popD_ind2                	popB_ind1                	0.2176982159169541       
DIST(1)          	476          	31                 	11                 	popD_ind2                	popB_ind2                	0.2259237146110073       
DIST(1)          	477          	31                 	12                 	popD_ind2                	popB_ind3                	0.27272862947566046      
DIST(1)          	478          	31                 	13                 	popD_ind2                	popB_ind4                	0.20245036636524588      
DIST(1)          	479          	31                 	14                 	popD_ind2                	popB_ind5                	0.25020944057727351      
DIST(1)          	480          	31                 	15                 	popD_ind2                	popB_ind6                	0.26265371451431879      
DIST(1)          	481          	31                 	16                 	popD_ind2                	popB_ind7                	0.22653330014344872      
DIST(1)          	482          	31                 	17                 	popD_ind2                	popB_ind8                	0.23948524980857966      
DIST(1)          	483          	31                 	18                 	popD_ind2                	popB_ind9                	0.21237658573738638      
DIST(1)          	484          	31                 	19                 	popD_ind2                	popB_ind10               	0.19899818404453432      
DIST(1)          	485          	31                 	20                 	popD_ind2                	popC_ind1                	0.16846483892402281      
DIST(1)          	486          	31                 	21                 	popD_ind2                	popC_ind2                	0.16672873162384189      
DIST(1)          	487          	31                 	22                 	popD_ind2                	popC_ind3                	0.16175438219767996      
DIST(1)          	488          	31                 	23                 	popD_ind2                	popC_ind4                	0.082048556450437946     
DIST(1)          	489          	31                 	24                 	popD_ind2                	popC_ind5                	0.065738464466928437     
DIST(1)          	490          	31                 	25                 	popD_ind2                	popC_ind6                	0.12302478124949576      
DIST(1)          	491          	31                 	26                 	popD_ind2                	popC_ind7                	0.16241923888461929      
DIST(1)          	492          	31                 	27                 	popD_ind2                	popC_ind8                	0.17524499344050826      
DIST(1)          	493          	31                 	28                 	popD_ind2                	popC_ind9                	0.11883894416167082      
DIST(1)          	494          	31                 	29                 	popD_ind2                	popC_ind10               	0.17443339134304414      
DIST(1)          	495          	31                 	30                 	popD_ind2                	popD_ind1                	0.064576161974153209     
DIST(1)          	496          	32                 	0                  	popD_ind3                	popA_ind1                	0.22101561676569043      
DIST(1)          	497          	32                 	1                  	popD_ind3                	popA_ind2                	0.091315068731080556     
DIST(1)          	498          	32                 	2                  	popD_ind3                	popA_ind3                	0.17047831749440687      
DIST(1)          	499          	32                 	3                  	popD_ind3                	popA_ind4                	0.18281749583573978      
DIST(1)          	500          	32                 	4                  	popD_ind3                	popA_ind5                	0.19268107311000743      
DIST(1)          	501          	32                 	5                  	popD_ind3                	popA_ind6                	0.22939261178266604      
DIST(1)          	502          	32                 	6                  	popD_ind3                	popA_ind7                	0.17338542867880868      
DIST(1)          	503          	32                 	7                  	popD_ind3                	popA_ind8                	0.19837001276276711      
DIST(1)          	504          	32                 	8                  	popD_ind3                	popA_ind9                	0.15008270114662803      
DIST(1)          	505          	32                 	9                  	popD_ind3                	popA_ind10               	0.18563201568799231      
DIST(1)          	506          	32                 	10                 	popD_ind3                	popB_ind1                	0.21893510830319898      
DIST(1)          	507          	32                 	11                 	popD_ind3                	popB_ind2                	0.18349102271350137      
DIST(1)          	508          	32                 	12                 	popD_ind3                	popB_ind3                	0.21277993095707076      
DIST(1)          	509          	32                 	13                 	popD_ind3                	popB_ind4                	0.17493860817319776      
DIST(1)          	510          	32                 	14                 	popD_ind3                	popB_ind5                	0.19796806509855139      
DIST(1)          	511          	32                 	15                 	popD_ind3                	popB_ind6                	0.17308135954369902      
DIST(1)          	512          	32                 	16                 	popD_ind3                	popB_ind7                	0.20618120362975101      
DIST(1)          	513          	32                 	17                 	popD_ind3                	popB_ind8                	0.14634320153739722      
DIST(1)          	514          	32                 	18                 	popD_ind3                	popB_ind9                	0.17480534081730836      
DIST(1)          	515          	32                 	19                 	popD_ind3                	popB_ind10               	0.17374194623019329      
DIST(1)          	516          	32                 	20                 	popD_ind3                	popC_ind1                	0.10930288882255607      
DIST(1)          	517          	32                 	21                 	popD_ind3                	popC_ind2                	0.10846395089436806      
DIST(1)          	518          	32                 	22                 	popD_ind3                	popC_ind3                	0.12783401100466402      
DIST(1)          	519          	32                 	23                 	popD_ind3                	popC_ind4                	0.114693748356072        
DIST(1)          	520          	32                 	24                 	popD_ind3                	popC_ind5                	0.08196968133652266      
DIST(1)          	521          	32                 	25                 	popD_ind3                	popC_ind6                	0.13636551983836262      
DIST(1)          	522          	32                 	26                 	popD_ind3                	popC_ind7                	0.16349176700277326      
DIST(1)          	523          	32                 	27                 	popD_ind3                	popC_ind8                	0.1192772196144207       
DIST(1)          	524          	32                 	28                 	popD_ind3                	popC_ind9                	0.12640741666249511      
DIST(1)          	525          	32                 	29                 	popD_ind3                	popC_ind10               	0.20027109942585403      
DIST(1)          	526          	32                 	30                 	popD_ind3                	popD_ind1                	0.055549174616355029     
DIST(1)          	527          	32                 	31                 	popD_ind3                	popD_ind2                	0.084741591710787267     
DIST(1)          	528          	33                 	0                  	popD_ind4                	popA_ind1                	0.20867651590808911      
DIST(1)          	529          	33                 	1                  	popD_ind4                	popA_ind2                	0.15860767321324087      
DIST(1)          	530          	33                 	2                  	popD_ind4                	popA_ind3                	0.20642891361451224      
DIST(1)          	531          	33                 	3                  	popD_ind4                	popA_ind4                	0.19888747055935213      
DIST(1)          	532          	33                 	4                  	popD_ind4                	popA_ind5                	0.1839241015249557       
DIST(1)          	533          	33                 	5                  	popD_ind4                	popA_ind6                	0.22144594242960741      
DIST(1)          	534          	33                 	6                  	popD_ind4                	popA_ind7                	0.17902726726923535      
DIST(1)          	535          	33                 	7                  	popD_ind4                	popA_ind8                	0.20990936913047711      
DIST(1)          	536          	33                 	8                  	popD_ind4                	popA_ind9                	0.19507804468020873      
DIST(1)          	537          	33                 	9                  	popD_ind4                	popA_ind10               	0.18038684264531968      
DIST(1)          	538          	33                 	10                 	popD_ind4                	popB_ind1                	0.23919289155923118      
DIST(1)          	539          	33                 	11                 	popD_ind4                	popB_ind2                	0.15478947561292922      
DIST(1)          	540          	33                 	12                 	popD_ind4                	popB_ind3                	0.19941233218382939      
DIST(1)          	541          	33                 	13                 	popD_ind4                	popB_ind4                	0.2057376434010979       
DIST(1)          	542          	33                 	14                 	popD_ind4                	popB_ind5                	0.25391447973445563      
DIST(1)          	543          	33                 	15                 	popD_ind4                	popB_ind6                	0.23530740789842941      
DIST(1)          	544          	33                 	16                 	popD_ind4                	popB_ind7                	0.22514748421643288      
DIST(1)          	545          	33                 	17                 	popD_ind4                	popB_ind8                	0.21471598894825328      
DIST(1)          	546          	33                 	18                 	popD_ind4                	popB_ind9                	0.19174464449363671      
DIST(1)          	547          	33                 	19                 	popD_ind4                	popB_ind10               	0.21651897983173196      
DIST(1)          	548          	33                 	20                 	popD_ind4                	popC_ind1                	0.14684020519339988      
DIST(1)          	549          	33                 	21                 	popD_ind4                	popC_ind2                	0.12154685529474629      
DIST(1)          	550          	33                 	22                 	popD_ind4                	popC_ind3                	0.16840734321912118      
DIST(1)          	551          	33                 	23                 	popD_ind4                	popC_ind4                	0.13102835165191384      
DIST(1)          	552          	33                 	24                 	popD_ind4                	popC_ind5                	0.10528296756241991      
DIST(1)          	553          	33                 	25                 	popD_ind4                	popC_ind6                	0.09245518859460268      
DIST(1)          	554          	33                 	26                 	popD_ind4                	popC_ind7                	0.15531964686520144      
DIST(1)          	555          	33                 	27                 	popD_ind4                	popC_ind8                	0.17810137419091554      
DIST(1)          	556          	33                 	28                 	popD_ind4                	popC_ind9                	0.11642410183818178      
DIST(1)          	557          	33                 	29                 	popD_ind4                	popC_ind10               	0.16329440387638294      
DIST(1)          	558          	33                 	30                 	popD_ind4                	popD_ind1                	0.040291159308572425     
DIST(1)          	559          	33                 	31                 	popD_ind4                	popD_ind2                	0.071182101174123102     
DIST(1)          	560          	33                 	32                 	popD_ind4                	popD_ind3                	0.050850306031042632     
DIST(1)          	561          	34                 	0                  	popD_ind5                	popA_ind1                	0.22885010737159645      
DIST(1)          	562          	34                 	1                  	popD_ind5                	popA_ind2                	0.23106126738994104      
DIST(1)          	563          	34                 	2                  	popD_ind5                	popA_ind3                	0.21736566698621618      
DIST(1)          	564          	34                 	3                  	popD_ind5                	popA_ind4                	0.21130195979700109      
DIST(1)          	565          	34                 	4                  	popD_ind5                	popA_ind5                	0.25371929719602032      
DIST(1)          	566          	34                 	5                  	popD_ind5                	popA_ind6                	0.2515874490593108       
DIST(1)          	567          	34                 	6                  	popD_ind5                	popA_ind7                	0.23206171111174129      
DIST(1)          	568          	34                 	7                  	popD_ind5                	popA_ind8                	0.22433396965282582      
DIST(1)          	569          	34                 	8                  	popD_ind5                	popA_ind9                	0.20455029119506379      
DIST(1)          	570          	34                 	9                  	popD_ind5                	popA_ind10               	0.2252583763083077       
DIST(1)          	571          	34                 	10                 	popD_ind5                	popB_ind1                	0.22001406799952866      
DIST(1)          	572          	34                 	11                 	popD_ind5                	popB_ind2                	0.1796267202123627       
DIST(1)          	573          	34                 	12                 	popD_ind5                	popB_ind3                	0.1907071028387467       
DIST(1)          	574          	34                 	13                 	popD_ind5                	popB_ind4                	0.17960029739689207      
DIST(1)          	575          	34                 	14                 	popD_ind5                	popB_ind5                	0.24667554847277948      
DIST(1)          	576          	34                 	15                 	popD_ind5                	popB_ind6                	0.24720539059618946      
DIST(1)          	577          	34                 	16                 	popD_ind5                	popB_ind7                	0.22253210946423457      
DIST(1)          	578          	34                 	17                 	popD_ind5                	popB_ind8                	0.22024858207637071      
DIST(1)          	579          	34                 	18                 	popD_ind5                	popB_ind9                	0.22087491212067956      
DIST(1)          	580          	34                 	19                 	popD_ind5                	popB_ind10               	0.20637398862038292      
DIST(1)          	581          	34                 	20                 	popD_ind5                	popC_ind1                	0.1702861277583908       
DIST(1)          	582          	34                 	21                 	popD_ind5                	popC_ind2                	0.13632496110792794      
DIST(1)          	583          	34                 	22                 	popD_ind5                	popC_ind3                	0.12864368350949137      
DIST(1)          	584          	34                 	23                 	popD_ind5                	popC_ind4                	0.1023151130244546       
DIST(1)          	585          	34                 	24                 	popD_ind5                	popC_ind5                	0.1011195831002315       
DIST(1)          	586          	34                 	25                 	popD_ind5                	popC_ind6                	0.10671319779168148      
DIST(1)          	587          	34                 	26                 	popD_ind5                	popC_ind7                	0.098049834520655646     
DIST(1)          	588          	34                 	27                 	popD_ind5                	popC_ind8                	0.1649475930962705       
DIST(1)          	589          	34                 	28                 	popD_ind5                	popC_ind9                	0.12493239032732528      
DIST(1)          	590          	34                 	29                 	popD_ind5                	popC_ind10               	0.19108624117004619      
DIST(1)          	591          	34                 	30                 	popD_ind5                	popD_ind1                	0.10955855872080915      
DIST(1)          	592          	34                 	31                 	popD_ind5                	popD_ind2                	0.063762947474836018     
DIST(1)          	593          	34                 	32                 	popD_ind5                	popD_ind3                	0.09797383948466136      
DIST(1)          	594          	34                 	33                 	popD_ind5                	popD_ind4                	0.070864105330816121     
DIST(1)          	595          	35                 	0                  	popD_ind6                	popA_ind1                	0.18906861575224948      
DIST(1)          	596          	35                 	1                  	popD_ind6                	popA_ind2                	0.17080634713491763      
DIST(1)          	597          	35                 	2                  	popD_ind6                	popA_ind3                	0.15972852315450686      
DIST(1)          	598          	35                 	3                  	popD_ind6                	popA_ind4                	0.19317388248819745      
DIST(1)          	599          	35                 	4                  	popD_ind6                	popA_ind5                	0.21104250591067072      
DIST(1)          	600          	35                 	5                  	popD_ind6                	popA_ind6                	0.23324782242385483      
DIST(1)          	601          	35                 	6                  	popD_ind6                	popA_ind7                	0.21248928031724801      
DIST(1)          	602          	35                 	7                  	popD_ind6                	popA_ind8                	0.19155181086985462      
DIST(1)          	603          	35                 	8                  	popD_ind6                	popA_ind9                	0.17139140749281642      
DIST(1)          	604          	35                 	9                  	popD_ind6                	popA_ind10               	0.19303624643676875      
DIST(1)          	605          	35                 	10                 	popD_ind6                	popB_ind1                	0.17745668766928208      
DIST(1)          	606          	35                 	11                 	popD_ind6                	popB_ind2                	0.16902646191598483      
DIST(1)          	607          	35                 	12                 	popD_ind6                	popB_ind3                	0.22158322180145551      
DIST(1)          	608          	35                 	13                 	popD_ind6                	popB_ind4                	0.20521023961175655      
DIST(1)          	609          	35                 	14                 	popD_ind6                	popB_ind5                	0.24950114791249173      
DIST(1)          	610          	35                 	15                 	popD_ind6                	popB_ind6                	0.27602586879920377      
DIST(1)          	611          	35                 	16                 	popD_ind6                	popB_ind7                	0.17684917109568893      
DIST(1)          	612          	35                 	17                 	popD_ind6                	popB_ind8                	0.17799709249866907      
DIST(1)          	613          	35                 	18                 	popD_ind6                	popB_ind9                	0.21912849783013125      
DIST(1)          	614          	35                 	19                 	popD_ind6                	popB_ind10               	0.19227054846838759      
DIST(1)          	615          	35                 	20                 	popD_ind6                	popC_ind1                	0.16985087038204127      
DIST(1)          	616          	35                 	21                 	popD_ind6                	popC_ind2                	0.13574686785218687      
DIST(1)          	617          	35                 	22                 	popD_ind6                	popC_ind3                	0.14363522563370673      
DIST(1)          	618          	35                 	23                 	popD_ind6                	popC_ind4                	0.14073379222291149      
DIST(1)          	619          	35                 	24                 	popD_ind6                	popC_ind5                	0.15605931194420619      
DIST(1)          	620          	35                 	25                 	popD_ind6                	popC_ind6                	0.13155786527481306      
DIST(1)          	621          	35                 	26                 	popD_ind6                	popC_ind7                	0.1376992537806139       
DIST(1)          	622          	35                 	27                 	popD_ind6                	popC_ind8                	0.17670211034180558      
DIST(1)          	623          	35                 	28                 	popD_ind6                	popC_ind9                	0.15302711426656274      
DIST(1)          	624          	35                 	29                 	popD_ind6                	popC_ind10               	0.20836779159024724      
DIST(1)          	625          	35                 	30                 	popD_ind6                	popD_ind1                	0.058841820962814241     
DIST(1)          	626          	35                 	31                 	popD_ind6                	popD_ind2                	0.085924998986668777     
DIST(1)          	627          	35                 	32                 	popD_ind6                	popD_ind3                	0.055913990614503112     
DIST(1)          	628          	35                 	33                 	popD_ind6                	popD_ind4                	0.054278446686643847     
DIST(1)          	629          	35                 	34                 	popD_ind6                	popD_ind5                	0.08341954273659806      
DIST(1)          	630          	36                 	0                  	popD_ind7                	popA_ind1                	0.21748781832871283      
DIST(1)          	631          	36                 	1                  	popD_ind7                	popA_ind2                	0.18133057044880249      
DIST(1)          	632          	36                 	2                  	popD_ind7                	popA_ind3                	0.15182524495064084      
DIST(1)          	633          	36                 	3                  	popD_ind7                	popA_ind4                	0.22336096806823813      
DIST(1)          	634          	36                 	4                  	popD_ind7                	popA_ind5                	0.23265503843116542      
DIST(1)          	635          	36                 	5                  	popD_ind7                	popA_ind6                	0.204857776332975        
DIST(1)          	636          	36                 	6                  	popD_ind7                	popA_ind7                	0.24349038518782767      
DIST(1)          	637          	36                 	7                  	popD_ind7                	popA_ind8                	0.18296823283842908      
DIST(1)          	638          	36                 	8                  	popD_ind7                	popA_ind9                	0.15268615997602997      
DIST(1)          	639          	36                 	9                  	popD_ind7                	popA_ind10               	0.2414335536970677       
DIST(1)          	640          	36                 	10                 	popD_ind7                	popB_ind1                	0.2310193038070813       
DIST(1)          	641          	36                 	11                 	popD_ind7                	popB_ind2                	0.12728144468959007      
DIST(1)          	642          	36                 	12                 	popD_ind7                	popB_ind3                	0.15930669423096996      
DIST(1)          	643          	36                 	13                 	popD_ind7                	popB_ind4                	0.22072638825434257      
DIST(1)          	644          	36                 	14                 	popD_ind7                	popB_ind5                	0.25550578838665883      
DIST(1)          	645          	36                 	15                 	popD_ind7                	popB_ind6                	0.20785729338976966      
DIST(1)          	646          	36                 	16                 	popD_ind7                	popB_ind7                	0.23757443620276592      
DIST(1)          	647          	36                 	17                 	popD_ind7                	popB_ind8                	0.21913313280021357      
DIST(1)          	648          	36                 	18                 	popD_ind7                	popB_ind9                	0.24753815852148192      
DIST(1)          	649          	36                 	19                 	popD_ind7                	popB_ind10               	0.20975835786253622      
DIST(1)          	650          	36                 	20                 	popD_ind7                	popC_ind1                	0.1612283372648694       
DIST(1)          	651          	36                 	21                 	popD_ind7                	popC_ind2                	0.13455722684136431      
DIST(1)          	652          	36                 	22                 	popD_ind7                	popC_ind3                	0.15821334044127192      
DIST(1)          	653          	36                 	23                 	popD_ind7                	popC_ind4                	0.17025141661988222      
DIST(1)          	654          	36                 	24                 	popD_ind7                	popC_ind5                	0.13967509150158611      
DIST(1)          	655          	36                 	25                 	popD_ind7                	popC_ind6                	0.17441079092696882      
DIST(1)          	656          	36                 	26                 	popD_ind7                	popC_ind7                	0.16216400083428653      
DIST(1)          	657          	36                 	27                 	popD_ind7                	popC_ind8                	0.17065790589462426      
DIST(1)          	658          	36                 	28                 	popD_ind7                	popC_ind9                	0.19962883261429851      
DIST(1)          	659          	36                 	29                 	popD_ind7                	popC_ind10               	0.18540859914235114      
DIST(1)          	660          	36                 	30                 	popD_ind7                	popD_ind1                	0.081386439574588451     
DIST(1)          	661          	36                 	31                 	popD_ind7                	popD_ind2                	0.05928901619108113      
DIST(1)          	662          	36                 	32                 	popD_ind7                	popD_ind3                	0.066235445652526467     
DIST(1)          	663          	36                 	33                 	popD_ind7                	popD_ind4                	0.052475010496292794     
DIST(1)          	664          	36                 	34                 	popD_ind7                	popD_ind5                	0.10335837115998307      
DIST(1)          	665          	36                 	35                 	popD_ind7                	popD_ind6                	0.10181903563687003      
DIST(1)          	666          	37                 	0                  	popD_ind8                	popA_ind1                	0.28611408867877525      
DIST(1)          	667          	37                 	1                  	popD_ind8                	popA_ind2                	0.19701157640374301      
DIST(1)          	668          	37                 	2                  	popD_ind8                	popA_ind3                	0.20851366427915058      
DIST(1)          	669          	37                 	3                  	popD_ind8                	popA_ind4                	0.2175316212354968       
DIST(1)          	670          	37                 	4                  	popD_ind8                	popA_ind5                	0.23075670961912625      
DIST(1)          	671          	37                 	5                  	popD_ind8                	popA_ind6                	0.26979832917824487      
DIST(1)          	672          	37                 	6                  	popD_ind8                	popA_ind7                	0.24279219748868716      
DIST(1)          	673          	37                 	7                  	popD_ind8                	popA_ind8                	0.18677507061986898      
DIST(1)          	674          	37                 	8                  	popD_ind8                	popA_ind9                	0.23162921250311158      
DIST(1)          	675          	37                 	9                  	popD_ind8                	popA_ind10               	0.26385518040665795      
DIST(1)          	676          	37                 	10                 	popD_ind8                	popB_ind1                	0.2318025252626327       
DIST(1)          	677          	37                 	11                 	popD_ind8                	popB_ind2                	0.25646514661602066      
DIST(1)          	678          	37                 	12                 	popD_ind8                	popB_ind3                	0.22949421606864762      
DIST(1)          	679          	37                 	13                 	popD_ind8                	popB_ind4                	0.21021289975470173      
DIST(1)          	680          	37                 	14                 	popD_ind8                	popB_ind5                	0.24035674846257959      
DIST(1)          	681          	37                 	15                 	popD_ind8                	popB_ind6                	0.23246961300340263      
DIST(1)          	682          	37                 	16                 	popD_ind8                	popB_ind7                	0.24919367820131316      
DIST(1)          	683          	37                 	17                 	popD_ind8                	popB_ind8                	0.23188610409305624      
DIST(1)          	684          	37                 	18                 	popD_ind8                	popB_ind9                	0.21870880054468225      
DIST(1)          	685          	37                 	19                 	popD_ind8                	popB_ind10               	0.2750491602878905       
DIST(1)          	686          	37                 	20                 	popD_ind8                	popC_ind1                	0.14123894476618126      
DIST(1)          	687          	37                 	21                 	popD_ind8                	popC_ind2                	0.14543631086115855      
DIST(1)          	688          	37                 	22                 	popD_ind8                	popC_ind3                	0.14769445662723019      
DIST(1)          	689          	37                 	23                 	popD_ind8                	popC_ind4                	0.1324365005803228       
DIST(1)          	690          	37                 	24                 	popD_ind8                	popC_ind5                	0.13959193467527084      
DIST(1)          	691          	37                 	25                 	popD_ind8                	popC_ind6                	0.18151937177781802      
DIST(1)          	692          	37                 	26                 	popD_ind8                	popC_ind7                	0.21248647886909616      
DIST(1)          	693          	37                 	27                 	popD_ind8                	popC_ind8                	0.13484168233253779      
DIST(1)          	694          	37                 	28                 	popD_ind8                	popC_ind9                	0.15687661938295838      
DIST(1)          	695          	37                 	29                 	popD_ind8                	popC_ind10               	0.1245720895951411       
DIST(1)          	696          	37                 	30                 	popD_ind8                	popD_ind1                	0.078822163626686245     
DIST(1)          	697          	37                 	31                 	popD_ind8                	popD_ind2                	0.031588077619130926     
DIST(1)          	698          	37                 	32                 	popD_ind8                	popD_ind3                	0.11223839620790221      
DIST(1)          	699          	37                 	33                 	popD_ind8                	popD_ind4                	0.1012011966232126       
DIST(1)          	700          	37                 	34                 	popD_ind8                	popD_ind5                	0.11154495870553248      
DIST(1)          	701          	37                 	35                 	popD_ind8                	popD_ind6                	0.12718764441321861      
DIST(1)          	702          	37                 	36                 	popD_ind8                	popD_ind7                	0.080805289358370766     
DIST(1)          	703          	38                 	0                  	popD_ind9                	popA_ind1                	0.19386377200036509      
DIST(1)          	704          	38                 	1                  	popD_ind9                	popA_ind2                	0.17041531683237979      
DIST(1)          	705          	38                 	2                  	popD_ind9                	popA_ind3                	0.18011104860037708      
DIST(1)          	706          	38                 	3                  	popD_ind9                	popA_ind4                	0.18727158297348423      
DIST(1)          	707          	38                 	4                  	popD_ind9                	popA_ind5                	0.17363076858038587      
DIST(1)          	708          	38                 	5                  	popD_ind9                	popA_ind6                	0.18137935721932286      
DIST(1)          	709          	38                 	6                  	popD_ind9                	popA_ind7                	0.15964983626200091      
DIST(1)          	710          	38                 	7                  	popD_ind9                	popA_ind8                	0.19857084662885979      
DIST(1)          	711          	38                 	8                  	popD_ind9                	popA_ind9                	0.1896546613283232       
DIST(1)          	712          	38                 	9                  	popD_ind9                	popA_ind10               	0.18684575168150394      
DIST(1)          	713          	38                 	10                 	popD_ind9                	popB_ind1                	0.21999460631057782      
DIST(1)          	714          	38                 	11                 	popD_ind9                	popB_ind2                	0.16089735821459797      
DIST(1)          	715          	38                 	12                 	popD_ind9                	popB_ind3                	0.1618335652338414       
DIST(1)          	716          	38                 	13                 	popD_ind9                	popB_ind4                	0.17926301699703348      
DIST(1)          	717          	38                 	14                 	popD_ind9                	popB_ind5                	0.20773487075049937      
DIST(1)          	718          	38                 	15                 	popD_ind9                	popB_ind6                	0.16792936705061307      
DIST(1)          	719          	38                 	16                 	popD_ind9                	popB_ind7                	0.19351853098772492      
DIST(1)          	720          	38                 	17                 	popD_ind9                	popB_ind8                	0.1978197016431148       
DIST(1)          	721          	38                 	18                 	popD_ind9                	popB_ind9                	0.20858438752073355      
DIST(1)          	722          	38                 	19                 	popD_ind9                	popB_ind10               	0.24375939834005073      
DIST(1)          	723          	38                 	20                 	popD_ind9                	popC_ind1                	0.14987999185384393      
DIST(1)          	724          	38                 	21                 	popD_ind9                	popC_ind2                	0.10901374432887036      
DIST(1)          	725          	38                 	22                 	popD_ind9                	popC_ind3                	0.17818709553400208      
DIST(1)          	726          	38                 	23                 	popD_ind9                	popC_ind4                	0.11457855711561778      
DIST(1)          	727          	38                 	24                 	popD_ind9                	popC_ind5                	0.10746869968649357      
DIST(1)          	728          	38                 	25                 	popD_ind9                	popC_ind6                	0.1353637847549248       
DIST(1)          	729          	38                 	26                 	popD_ind9                	popC_ind7                	0.19869567622299128      
DIST(1)          	730          	38                 	27                 	popD_ind9                	popC_ind8                	0.12890692492276618      
DIST(1)          	731          	38                 	28                 	popD_ind9                	popC_ind9                	0.10481242462390886      
DIST(1)          	732          	38                 	29                 	popD_ind9                	popC_ind10               	0.11043098197863571      
DIST(1)          	733          	38                 	30                 	popD_ind9                	popD_ind1                	0.085772701730966169     
DIST(1)          	734          	38                 	31                 	popD_ind9                	popD_ind2                	0.061212496725935261     
DIST(1)          	735          	38                 	32                 	popD_ind9                	popD_ind3                	0.10084863335972481      
DIST(1)          	736          	38                 	33                 	popD_ind9                	popD_ind4                	0.044549299427491008     
DIST(1)          	737          	38                 	34                 	popD_ind9                	popD_ind5                	0.13737949274062225      
DIST(1)          	738          	38                 	35                 	popD_ind9                	popD_ind6                	0.117635130147692        
DIST(1)          	739          	38                 	36                 	popD_ind9                	popD_ind7                	0.07376582539050959      
DIST(1)          	740          	38                 	37                 	popD_ind9                	popD_ind8                	0.058150291609240699     
DIST(1)          	741          	39                 	0                  	popD_ind10               	popA_ind1                	0.23591998461270303      
DIST(1)          	742          	39                 	1                  	popD_ind10               	popA_ind2                	0.18141088869471156      
DIST(1)          	743          	39                 	2                  	popD_ind10               	popA_ind3                	0.17512422718603643      
DIST(1)          	744          	39                 	3                  	popD_ind10               	popA_ind4                	0.21880728377739112      
DIST(1)          	745          	39                 	4                  	popD_ind10               	popA_ind5                	0.23389374940547469      
DIST(1)          	746          	39                 	5                  	popD_ind10               	popA_ind6                	0.25865775635084864      
DIST(1)          	747          	39                 	6                  	popD_ind10               	popA_ind7                	0.25736482074151124      
DIST(1)          	748          	39                 	7                  	popD_ind10               	popA_ind8                	0.2258648684413308       
DIST(1)          	749          	39                 	8                  	popD_ind10               	popA_ind9                	0.19813393531297313      
DIST(1)          	750          	39                 	9                  	popD_ind10               	popA_ind10               	0.24799146793248933      
DIST(1)          	751          	39                 	10                 	popD_ind10               	popB_ind1                	0.2695020158895079       
DIST(1)          	752          	39                 	11                 	popD_ind10               	popB_ind2                	0.18603445591209467      
DIST(1)          	753          	39                 	12                 	popD_ind10               	popB_ind3                	0.22598027326955608      
DIST(1)          	754          	39                 	13                 	popD_ind10               	popB_ind4                	0.24420435505452867      
DIST(1)          	755          	39                 	14                 	popD_ind10               	popB_ind5                	0.26915977775154021      
DIST(1)          	756          	39                 	15                 	popD_ind10               	popB_ind6                	0.26054458549968329      
DIST(1)          	757          	39                 	16                 	popD_ind10               	popB_ind7                	0.23159659610687197      
DIST(1)          	758          	39                 	17                 	popD_ind10               	popB_ind8                	0.25040070451259233      
DIST(1)          	759          	39                 	18                 	popD_ind10               	popB_ind9                	0.22337389092966403      
DIST(1)          	760          	39                 	19                 	popD_ind10               	popB_ind10               	0.20537965524160698      
DIST(1)          	761          	39                 	20                 	popD_ind10               	popC_ind1                	0.18979746442769108      
DIST(1)          	762          	39                 	21                 	popD_ind10               	popC_ind2                	0.1432743210487315       
DIST(1)          	763          	39                 	22                 	popD_ind10               	popC_ind3                	0.19537235831398669      
DIST(1)          	764          	39                 	23                 	popD_ind10               	popC_ind4                	0.14155429190999963      
DIST(1)          	765          	39                 	24                 	popD_ind10               	popC_ind5                	0.12936714956764023      
DIST(1)          	766          	39                 	25                 	popD_ind10               	popC_ind6                	0.13009028477779366      
DIST(1)          	767          	39                 	26                 	popD_ind10               	popC_ind7                	0.15849379761903853      
DIST(1)          	768          	39                 	27                 	popD_ind10               	popC_ind8                	0.1362987691971379       
DIST(1)          	769          	39                 	28                 	popD_ind10               	popC_ind9                	0.15005500084582202      
DIST(1)          	770          	39                 	29                 	popD_ind10               	popC_ind10               	0.19704344649624525      
DIST(1)          	771          	39                 	30                 	popD_ind10               	popD_ind1                	0.095958891846744526     
DIST(1)          	772          	39                 	31                 	popD_ind10               	popD_ind2                	0.076421736302500695     
DIST(1)          	773          	39                 	32                 	popD_ind10               	popD_ind3                	0.076350578502792588     
DIST(1)          	774          	39                 	33                 	popD_ind10               	popD_ind4                	0.1039593327513201       
DIST(1)          	775          	39                 	34                 	popD_ind10               	popD_ind5                	0.075444102704649435     
DIST(1)          	776          	39                 	35                 	popD_ind10               	popD_ind6                	0.097702750023059101     
DIST(1)          	777          	39                 	36                 	popD_ind10               	popD_ind7                	0.093001479590568689     
DIST(1)          	778          	39                 	37                 	popD_ind10               	popD_ind8                	0.070368147500554423     
DIST(1)          	779          	39                 	38                 	popD_ind10               	popD_ind9                	0.12730920142388133      
# DMAT(MATRIX_INDEX), Distance matrix in matrix format:
# DMAT(1)        	[2]ITEMS                 	popA_ind1                	popA_ind2                	popA_ind3                	popA_ind4                	popA_ind5                	popA_ind6                	popA_ind7                	popA_ind8                	popA_ind9                	popA_ind10               	popB_ind1                	popB_ind2                	popB_ind3                	popB_ind4                	popB_ind5                	popB_ind6                	popB_ind7                	popB_ind8                	popB_ind9                	popB_ind10               	popC_ind1                	popC_ind2                	popC_ind3                	popC_ind4                	popC_ind5                	popC_ind6                	popC_ind7                	popC_ind8                	popC_ind9                	popC_ind10               	popD_ind1                	popD_ind2                	popD_ind3                	popD_ind4                	popD_ind5                	popD_ind6                	popD_ind7                	popD_ind8                	popD_ind9                	popD_ind10               
DMAT(1)          	popA_ind1                	
DMAT(1)          	popA_ind2                	0.11467617878074976      
DMAT(1)          	popA_ind3                	0.093092202947423805     	0.043429845634977267     
DMAT(1)          	popA_ind4                	0.099068543008061838     	0.14476271136798247      	0.092509288625741343     
DMAT(1)          	popA_ind5                	0.11548749475140628      	0.10235826812419882      	0.09521987109482824      	0.10924768687545468      
DMAT(1)          	popA_ind6                	0.10717339691795183      	0.08124585764909141      	0.10665280767883602      	0.099225375133700472     	0.084457781887820299     
DMAT(1)          	popA_ind7                	0.10530027483526966      	0.10952915292971828      	0.056029316909226939     	0.11165936466273696      	0.11752068019191368      	0.099017757040388121     
DMAT(1)          	popA_ind8                	0.097166560252935116     	0.087725183305923007     	0.035463864283368333     	0.10971550268283756      	0.052211145413247562     	0.080090237697028083     	0.081349577339648454     
DMAT(1)          	popA_ind9                	0.090441910531874198     	0.095867900816094934     	0.058549180362201331     	0.1026648719420763       	0.069954678492391426     	0.092754438818220694     	0.067255283583122355     	0.081902410789416333     
DMAT(1)          	popA_ind10               	0.12644992674147423      	0.0966360124762633       	0.10716580364420696      	0.101843573924255        	0.10530742687306807      	0.079663031092544578     	0.13024246593697633      	0.086512330866314727     	0.091715136836088729     
DMAT(1)          	popB_ind1                	0.2336180839991483       	0.18954759077640648      	0.20283944424946287      	0.23700842812959011      	0.21720537528782002      	0.16646532629003263      	0.17940566801556734      	0.18782183744554465      	0.17880978453375065      	0.1331576985459951       
DMAT(1)          	popB_ind2                	0.1557493052135199       	0.12902711380341356      	0.18546996912365249      	0.20095151449527437      	0.13921199302111537      	0.11863830608516811      	0.21865482444832432      	0.16953188761612992      	0.13355896813518423      	0.13685061378832802      	0.082160885097566277     
DMAT(1)          	popB_ind3                	0.18435734223985234      	0.18073172098884416      	0.12921894115135346      	0.1228474585695065       	0.18104606093364703      	0.19471725344017074      	0.1936455547548524       	0.11517994569806499      	0.1348604073845501       	0.19398220480765643      	0.16274535265075502      	0.13373797837552709      
DMAT(1)          	popB_ind4                	0.13365602036443922      	0.17171434159611618      	0.14385416636385351      	0.1362039295122727       	0.13844854551032734      	0.20507900555218048      	0.16674800759215763      	0.19158365508898345      	0.14481723658856768      	0.17216895148979977      	0.10146661607076184      	0.12618862661658065      	0.058739516198463106     
DMAT(1)          	popB_ind5                	0.21820072728406747      	0.14650675561445772      	0.19196571222817471      	0.17618931611687783      	0.16128517489385924      	0.15657767928982225      	0.17567679918954107      	0.16873546816436644      	0.19334619577350504      	0.14622659497598672      	0.076306207154471964     	0.066529877811036503     	0.14094098427547055      	0.096059387487813036     
DMAT(1)          	popB_ind6                	0.1384893437417391       	0.15612309594158832      	0.17415526086063793      	0.19717244117112803      	0.2018992075184749       	0.21586937475895562      	0.18762449822932747      	0.19803684293859389      	0.17157542180045365      	0.15205420861488694      	0.13700903211250814      	0.08776638529312078      	0.14859284705566389      	0.13156715752246448      	0.10742285964688825      
DMAT(1)          	popB_ind7                	0.17973240429294551      	0.14078806644125305      	0.091761684712605388     	0.14915398801606841      	0.16024652644638937      	0.16262642621955042      	0.1238278896798205       	0.12664600605320162      	0.11738137106945715      	0.15210423771205944      	0.082060786862680424     	0.13608893161450789      	0.080284930300773313     	0.069619761483182191     	0.13062814343793538      	0.14559357684436566      
DMAT(1)          	popB_ind8                	0.16567511205365207      	0.12640773566882468      	0.17088661897786564      	0.23133671988418381      	0.19225338752473931      	0.19726913409919913      	0.14794873930210076      	0.17227982719132295      	0.14959957714123034      	0.15249983598143488      	0.096185128202881956     	0.12356963919485967      	0.097396995678103074     	0.11707025754164227      	0.077202720953757867     	0.12856793997731336      	0.12011526516704352      
DMAT(1)          	popB_ind9                	0.12538366546380877      	0.15367870315230986      	0.10515210094250806      	0.15781540012958445      	0.12211001480676478      	0.13742681971106563      	0.11147969497017537      	0.14821969503273408      	0.095865949319817359     	0.15796362529584293      	0.11159237196376755      	0.10566166338278882      	0.11834799348704117      	0.09777857849635356      	0.092115999449317851     	0.11123602906816323      	0.079927548958855135     	0.097300091716630593     
DMAT(1)          	popB_ind10               	0.21307666193115798      	0.1703770676526504       	0.14639700483239174      	0.19214745347349027      	0.13061699069084554      	0.15657229526571761      	0.16985562811883512      	0.18171893336639938      	0.1671966477195756       	0.16347371290641377      	0.12687343944407431      	0.14430340760295579      	0.071379753531828252     	0.074913241182063356     	0.11966181426920416      	0.14635490020317699      	0.093136998894170472     	0.079366584664830125     	0.11427538896625043      
DMAT(1)          	popC_ind1                	0.2015332824908172       	0.21319717169203276      	0.19935098429448844      	0.22274743747405051      	0.13482631394144035      	0.22349667533920595      	0.21226372044518796      	0.15416149791707118      	0.14624866514036183      	0.19273603119366223      	0.18992584418755426      	0.13916576111370926      	0.21871802943394594      	0.23894372190687277      	0.17615682862512969      	0.21229738138034671      	0.25276437247717154      	0.22832392155407683      	0.20775044417262878      	0.14419441349830023      
DMAT(1)          	popC_ind2                	0.22010928878735497      	0.19596048314884862      	0.13348044864677666      	0.22768146183832752      	0.17210693161572238      	0.22974283767662584      	0.20868695076057764      	0.21980391670649738      	0.20436847500812061      	0.2377392532730476       	0.24110193894173057      	0.18274924278785168      	0.22482126011468134      	0.21303988719470013      	0.24704618511131377      	0.17241113608750086      	0.22743761668678641      	0.24500901223963265      	0.17401907242445566      	0.17641728721519209      	0.085929952914894581     
DMAT(1)          	popC_ind3                	0.20992072881488483      	0.19856025357945117      	0.18137186054747734      	0.22563648648046841      	0.19355061640016896      	0.20568966173104594      	0.18097502863792991      	0.19008316661633073      	0.20914506792198156      	0.19661696569243406      	0.20516226119894762      	0.17419158370645485      	0.223667855964443        	0.23600090782376576      	0.24198609467351978      	0.19888742525337749      	0.21567635282876937      	0.20019062125200154      	0.20930096488445007      	0.19961998820517379      	0.06458724038052234      	0.083556136378852577     
DMAT(1)          	popC_ind4                	0.18011973954276245      	0.17799342432036946      	0.14461749966416013      	0.19707754359541055      	0.081157311261563378     	0.19920225108790207      	0.23824256705359412      	0.12737489517049941      	0.13183984482082806      	0.16902390438926712      	0.16848696867103963      	0.1485571300127603       	0.22618623471090815      	0.19733227421567021      	0.21928494446825603      	0.14381946282788433      	0.18400919938333196      	0.2153539390760478       	0.17428643282842946      	0.1677931257830779       	0.088996394857584066     	0.0636691408407625       	0.087365595598285314     
DMAT(1)          	popC_ind5                	0.20325323130980855      	0.14702450642836418      	0.12378919856248129      	0.16490775905017879      	0.12756045244004607      	0.15275739577477193      	0.14767875278433951      	0.11134245564017298      	0.1545797465156945       	0.16103514198867441      	0.18153955539439667      	0.12645621735530049      	0.21219604150214644      	0.22434630118986323      	0.21254986108129204      	0.16881940483549152      	0.14723056087689204      	0.20735358341427781      	0.15812908243462728      	0.12767489028796097      	0.038373039060623483     	0.092059644589807504     	0.074413818107818774     	0.053726777929373526     
DMAT(1)          	popC_ind6                	0.20984031673621628      	0.17442468605602213      	0.25895711393091969      	0.2116313684226776       	0.16355731741209173      	0.26326882256131307      	0.26393090198139951      	0.168157258263574        	0.22714663644873589      	0.24610165008253934      	0.24285448097296566      	0.18794323710484434      	0.2442810931287237       	0.18087135708310473      	0.18653439343397099      	0.2058210967759192       	0.22320203943759701      	0.21083984955048851      	0.18504443247286231      	0.1909788577390174       	0.061336527554696978     	0.060911166662706975     	0.088665456548880514     	0.07408113418657028      	0.053431671951617882     
DMAT(1)          	popC_ind7                	0.23655857744765849      	0.24830740085268679      	0.22165859931328433      	0.18967820697965587      	0.20138998552086085      	0.28795560981523416      	0.26873541148604974      	0.19675588091388402      	0.22811461996768165      	0.20867271024605341      	0.24061109001033376      	0.19877923217784804      	0.21547008840665588      	0.19767769172169303      	0.27509638561703653      	0.25374463338126252      	0.17857892695742442      	0.22978234231581421      	0.22425103724399476      	0.19836430073578476      	0.080468945192271074     	0.077573278036358162     	0.096624243274473581     	0.078420141511245087     	0.085692851086862698     	0.094015717217159261     
DMAT(1)          	popC_ind8                	0.231717259781694        	0.22090372507239445      	0.2080977742654021       	0.19897184497235881      	0.20981466147567945      	0.23941607366271042      	0.20018477116636779      	0.17347278125565413      	0.20419846762542376      	0.21781060486206255      	0.22480874276002183      	0.27709420638695609      	0.22787476566118273      	0.19083714046180264      	0.30645922919112695      	0.23059273142635173      	0.23069574344760815      	0.26584523918957531      	0.22086612680901907      	0.199701679793234        	0.038367044363519689     	0.057603914036350297     	0.10002365391216403      	0.074811407190029849     	0.049541727827786659     	0.083376674419933502     	0.10231720437086897      
DMAT(1)          	popC_ind9                	0.27893889901517854      	0.21410309795138516      	0.19286296697013494      	0.18686749308825443      	0.15456903913146036      	0.21061187153154848      	0.22837868173009906      	0.16378460291296684      	0.17143383683045693      	0.14333689465982724      	0.24031571106900579      	0.19753342149877429      	0.24909517760238459      	0.21178624032341026      	0.23701875375463241      	0.20842888691417302      	0.27955586354207179      	0.29242445900054814      	0.25257668494528773      	0.195070428112468        	0.062541815271601547     	0.096198048711564496     	0.12801917742055624      	0.095481720882076709     	0.046700724870550984     	0.078023295539278517     	0.094692791475228286     	0.094952695531240541     
DMAT(1)          	popC_ind10               	0.26359692965102527      	0.24990291569905129      	0.24506414677891714      	0.27575413391200854      	0.24426070118934101      	0.22138250336297258      	0.22809973495139563      	0.22615354130733084      	0.24653119815044344      	0.22949346461278047      	0.20338124312427924      	0.23936921589388915      	0.24991241331889913      	0.22499187380586272      	0.27313820987458215      	0.25544194022353878      	0.27546492200799055      	0.26841693667885036      	0.21351601714660889      	0.19710585973442807      	0.084133763987661134     	0.089002603093127963     	0.090197021640572708     	0.070476984843123952     	0.069405838498509634     	0.081919128661087498     	0.10398249075452973      	0.088913855721502205     	0.061781873583093866     
DMAT(1)          	popD_ind1                	0.15576874376073427      	0.18354603047325474      	0.15548244825433163      	0.23636189166322452      	0.16065407838892959      	0.25146596833430707      	0.18426987979109441      	0.14242565452649741      	0.18239425024363826      	0.20453525675933965      	0.20477225865304255      	0.20888693041062803      	0.19351193017586943      	0.19991354311571102      	0.24126918206266984      	0.17394113001215425      	0.18805765501821284      	0.18541141583156168      	0.15971641445464069      	0.20158372573305069      	0.12984999798110181      	0.093924844356783097     	0.1495852717103146       	0.087817162380279229     	0.074551563154910003     	0.11993332403573874      	0.16517261220996163      	0.11194600643558583      	0.13796409469097354      	0.16068347267613489      
DMAT(1)          	popD_ind2                	0.2252919853265799       	0.15273144787650705      	0.23137037387074694      	0.21336392171997709      	0.18413506961360537      	0.20093618363776039      	0.17424317963019309      	0.22066812525444188      	0.16343172198789696      	0.26837875257604882      	0.2176982159169541       	0.2259237146110073       	0.27272862947566046      	0.20245036636524588      	0.25020944057727351      	0.26265371451431879      	0.22653330014344872      	0.23948524980857966      	0.21237658573738638      	0.19899818404453432      	0.16846483892402281      	0.16672873162384189      	0.16175438219767996      	0.082048556450437946     	0.065738464466928437     	0.12302478124949576      	0.16241923888461929      	0.17524499344050826      	0.11883894416167082      	0.17443339134304414      	0.064576161974153209     
DMAT(1)          	popD_ind3                	0.22101561676569043      	0.091315068731080556     	0.17047831749440687      	0.18281749583573978      	0.19268107311000743      	0.22939261178266604      	0.17338542867880868      	0.19837001276276711      	0.15008270114662803      	0.18563201568799231      	0.21893510830319898      	0.18349102271350137      	0.21277993095707076      	0.17493860817319776      	0.19796806509855139      	0.17308135954369902      	0.20618120362975101      	0.14634320153739722      	0.17480534081730836      	0.17374194623019329      	0.10930288882255607      	0.10846395089436806      	0.12783401100466402      	0.114693748356072        	0.08196968133652266      	0.13636551983836262      	0.16349176700277326      	0.1192772196144207       	0.12640741666249511      	0.20027109942585403      	0.055549174616355029     	0.084741591710787267     
DMAT(1)          	popD_ind4                	0.20867651590808911      	0.15860767321324087      	0.20642891361451224      	0.19888747055935213      	0.1839241015249557       	0.22144594242960741      	0.17902726726923535      	0.20990936913047711      	0.19507804468020873      	0.18038684264531968      	0.23919289155923118      	0.15478947561292922      	0.19941233218382939      	0.2057376434010979       	0.25391447973445563      	0.23530740789842941      	0.22514748421643288      	0.21471598894825328      	0.19174464449363671      	0.21651897983173196      	0.14684020519339988      	0.12154685529474629      	0.16840734321912118      	0.13102835165191384      	0.10528296756241991      	0.09245518859460268      	0.15531964686520144      	0.17810137419091554      	0.11642410183818178      	0.16329440387638294      	0.040291159308572425     	0.071182101174123102     	0.050850306031042632     
DMAT(1)          	popD_ind5                	0.22885010737159645      	0.23106126738994104      	0.21736566698621618      	0.21130195979700109      	0.25371929719602032      	0.2515874490593108       	0.23206171111174129      	0.22433396965282582      	0.20455029119506379      	0.2252583763083077       	0.22001406799952866      	0.1796267202123627       	0.1907071028387467       	0.17960029739689207      	0.24667554847277948      	0.24720539059618946      	0.22253210946423457      	0.22024858207637071      	0.22087491212067956      	0.20637398862038292      	0.1702861277583908       	0.13632496110792794      	0.12864368350949137      	0.1023151130244546       	0.1011195831002315       	0.10671319779168148      	0.098049834520655646     	0.1649475930962705       	0.12493239032732528      	0.19108624117004619      	0.10955855872080915      	0.063762947474836018     	0.09797383948466136      	0.070864105330816121     
DMAT(1)          	popD_ind6                	0.18906861575224948      	0.17080634713491763      	0.15972852315450686      	0.19317388248819745      	0.21104250591067072      	0.23324782242385483      	0.21248928031724801      	0.19155181086985462      	0.17139140749281642      	0.19303624643676875      	0.17745668766928208      	0.16902646191598483      	0.22158322180145551      	0.20521023961175655      	0.24950114791249173      	0.27602586879920377      	0.17684917109568893      	0.17799709249866907      	0.21912849783013125      	0.19227054846838759      	0.16985087038204127      	0.13574686785218687      	0.14363522563370673      	0.14073379222291149      	0.15605931194420619      	0.13155786527481306      	0.1376992537806139       	0.17670211034180558      	0.15302711426656274      	0.20836779159024724      	0.058841820962814241     	0.085924998986668777     	0.055913990614503112     	0.054278446686643847     	0.08341954273659806      
DMAT(1)          	popD_ind7                	0.21748781832871283      	0.18133057044880249      	0.15182524495064084      	0.22336096806823813      	0.23265503843116542      	0.204857776332975        	0.24349038518782767      	0.18296823283842908      	0.15268615997602997      	0.2414335536970677       	0.2310193038070813       	0.12728144468959007      	0.15930669423096996      	0.22072638825434257      	0.25550578838665883      	0.20785729338976966      	0.23757443620276592      	0.21913313280021357      	0.24753815852148192      	0.20975835786253622      	0.1612283372648694       	0.13455722684136431      	0.15821334044127192      	0.17025141661988222      	0.13967509150158611      	0.17441079092696882      	0.16216400083428653      	0.17065790589462426      	0.19962883261429851      	0.18540859914235114      	0.081386439574588451     	0.05928901619108113      	0.066235445652526467     	0.052475010496292794     	0.10335837115998307      	0.10181903563687003      
DMAT(1)          	popD_ind8                	0.28611408867877525      	0.19701157640374301      	0.20851366427915058      	0.2175316212354968       	0.23075670961912625      	0.26979832917824487      	0.24279219748868716      	0.18677507061986898      	0.23162921250311158      	0.26385518040665795      	0.2318025252626327       	0.25646514661602066      	0.22949421606864762      	0.21021289975470173      	0.24035674846257959      	0.23246961300340263      	0.24919367820131316      	0.23188610409305624      	0.21870880054468225      	0.2750491602878905       	0.14123894476618126      	0.14543631086115855      	0.14769445662723019      	0.1324365005803228       	0.13959193467527084      	0.18151937177781802      	0.21248647886909616      	0.13484168233253779      	0.15687661938295838      	0.1245720895951411       	0.078822163626686245     	0.031588077619130926     	0.11223839620790221      	0.1012011966232126       	0.11154495870553248      	0.12718764441321861      	0.080805289358370766     
DMAT(1)          	popD_ind9                	0.19386377200036509      	0.17041531683237979      	0.18011104860037708      	0.18727158297348423      	0.17363076858038587      	0.18137935721932286      	0.15964983626200091      	0.19857084662885979      	0.1896546613283232       	0.18684575168150394      	0.21999460631057782      	0.16089735821459797      	0.1618335652338414       	0.17926301699703348      	0.20773487075049937      	0.16792936705061307      	0.19351853098772492      	0.1978197016431148       	0.20858438752073355      	0.24375939834005073      	0.14987999185384393      	0.10901374432887036      	0.17818709553400208      	0.11457855711561778      	0.10746869968649357      	0.1353637847549248       	0.19869567622299128      	0.12890692492276618      	0.10481242462390886      	0.11043098197863571      	0.085772701730966169     	0.061212496725935261     	0.10084863335972481      	0.044549299427491008     	0.13737949274062225      	0.117635130147692        	0.07376582539050959      	0.058150291609240699     
DMAT(1)          	popD_ind10               	0.23591998461270303      	0.18141088869471156      	0.17512422718603643      	0.21880728377739112      	0.23389374940547469      	0.25865775635084864      	0.25736482074151124      	0.2258648684413308       	0.19813393531297313      	0.24799146793248933      	0.2695020158895079       	0.18603445591209467      	0.22598027326955608      	0.24420435505452867      	0.26915977775154021      	0.26054458549968329      	0.23159659610687197      	0.25040070451259233      	0.22337389092966403      	0.20537965524160698      	0.18979746442769108      	0.1432743210487315       	0.19537235831398669      	0.14155429190999963      	0.12936714956764023      	0.13009028477779366      	0.15849379761903853      	0.1362987691971379       	0.15005500084582202      	0.19704344649624525      	0.095958891846744526     	0.076421736302500695     	0.076350578502792588     	0.1039593327513201       	0.075444102704649435     	0.097702750023059101     	0.093001479590568689     	0.070368147500554423     	0.12730920142388133      
```

[!TIP] As `--print-dm` is an INT+ argument, using `--print-dm 3` (i.e. 1+2) will print the distance matrix in both condensed and verbose format.

## Distance matrix file input 

Distance matrix file input option `--in-dm <filename>` expects a condensed distance matrix file (produced by ngsAMOVA `--print-dm 1`).
