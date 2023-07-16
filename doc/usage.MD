# Main analysis commands


## Actions (`-doAction <int>`)
Use action commands to specify the action (i.e. analysis) to be performed.
Action commands are of the form `-doXXXX <int>`
where `XXXX`` defines the general type of analysis to be performed and <int> specifies the exact analysis to be performed.

e.g. `-doAMOVA 1` specifies to perform AMOVA analysis with the genotype likelihoods.

The following action commands are available:
- `doAMOVA <int>` : perform AMOVA analysis
- `doEM <int>` : perform EM optimization
- `doDxy <int>` : estimate Dxy
- `doPhylo <int>` : do neighbor-joining
- `doDist <int>` : estimate pairwise distance matrix


| Command | Description | Requires | Argument values | Output |
| ------- | ----------- | -------- | ------- | ------ |
| `-doAMOVA` | Perform Analysis of Molecular Variance (AMOVA) | - Distance Matrix (`-doDist` or `--in-dm` ) <br/> - Metadata file (`--metadata/-m <FILE>`) <br/>- Formula (`--formula/-f <string>`) | [0:default] Do not run AMOVA <br/> [1] Run AMOVA with genotype likelihoods <br/> [2] Run AMOVA with genotypes | `<output_prefix>.amova.csv` |
| `-doDxy` | Calculate Absolute Nucleotide Divergence ($d_{XY}$) | - Distance Matrix (`-doDist` or `--in-dm`) | [0:default] Do not run Dxy <br/> [1] Run Dxy for all groups in each hierarchical level defined in the metadata file <br/> [2] Run Dxy for the groups in metadata specified by the user via `--dxy-groups` <br/> [3] Run Dxy for the groups in the hierarchical level specified by the user via `--dxy-levels` | `<output_prefix>.dxy.csv` |
| `-doPhylo` | Construct Phylogenetic Tree | - Distance Matrix (`-doDist` or `--in-dm`) | [0:default] Do not run phylogenetic tree construction <br/> [1] Construct phylogenetic tree using Neighbor-Joining with individuals as leaf nodes <br/> [2] Construct phylogenetic tree using Neighbor-Joining with groups as leaf nodes (requires: `-doDxy`) | `<output_prefix>.newick` |
| `-doDist` | Estimate pairwise distance matrix | Input VCF file with Genotype Likelihoods (GL tag) or genotypes (GT tag) (`--in-vcf/-i`) | [0:default] Do not estimate distance matrix <br/> [1] Estimate distance matrix using genotype likelihoods (requires: `-doEM`) <br/> [2] Estimate distance matrix using genotypes  <br/> [3] Read distance matrix from file | `<output_prefix>.distance_matrix.csv` |

___

## Define input and output files

### Input files (`--in-XXX`)

Specify the input file types and their filenames.
Input file is defined as the file containing the data to be used in the analyses.
Input file commands are of the form `--in-XXX <filename>` where XXX defines the type of input file and `<filename>` specifies the name of the input file.
For example, `--in-vcf <filename>` specifies to read the input file as a VCF file.

The following input file commands are available:
  - `--in-vcf <filename>` : VCF file input
  - `--in-dm <filename>` : Distance matrix input

### Output files (`--output/-o`)

Specify the output filename prefix.
Output file names are automatically generated based on the value of the `--output` argument and the file type to be printed.

For example, if `--output` is set to `out`, then the output file name for the AMOVA results will be `out.amova.csv`.

___

## Arguments (`--long-form/-sf`)

Use argument commands to specify the parameters to be used in the analyses.

Argument commands are of the form `--long-form <int>` or `-sf <int>`
  where long form of the argument starts with double dash `--` and separated by hyphen `-`
  and short form of the argument starts with single dash `-` and is typically the first letter(s) of the long form

| Command | Description | Requires | Argument values | Output |
| ------- | ----------- | -------- | ------- | ------ |
| `--dxy-groups` | Specify the groups to be used in the Dxy analysis | - Metadata file (`--metadata/-m <FILE>`) | Comma-separated list of group names | `<output_prefix>.dxy.csv` |
| `--dxy-levels` | Specify the hierarchical levels to be used in the Dxy analysis | - Metadata file (`--metadata/-m <FILE>`) | Comma-separated list of hierarchical levels | `<output_prefix>.dxy.csv` |
| `--seed` | Specify the seed for the random number generator | - | Integer | - |
| `--metadata` | Specify the metadata file | - | Metadata file | - |
| `--verbose` | Specify the verbosity level | - | Integer | - |
| `--formula` | Specify the formula for the AMOVA analysis | - | Formula | - |
| `--majorMinorFile` | Specify the major and minor alleles file | - | Major and minor alleles file | - |
| `--ancDerFile` | Specify the ancestral and derived alleles file | - | Ancestral and derived alleles file | - |

___

### Printing commands (`--printXxXxx/--printXX/-pXX`)

Use printing commands to specify the output files to be generated.

Printing commands are of the form `--printXxXxx <int>` or `--printXX <int>` or `-pXX <int>` where `XxXxx` is the long form of the file type to be printed, `XX` is the short form of the file type (typically  the first letter(s) of the long form), and `<int>` defines the compression level of the output file.

Output file names are automatically generated based on the value of the `--output` argument and the file type to be printed.

The following printing commands are available:
  - `--printJointGenotypeCountMatrix/-pJGCD <int>` : print joint genotype count distribution
  - `--printAmovaTable/-pAT <int>` : print AMOVA table
  - `--printDistanceMatrix/-pDM <int>` : print distance matrix
  - `--printBlocksTab <0|1>` : print tab-delimited blocks file defining the start and end positions of each block (default: 0 = do not print, 1 = print)
