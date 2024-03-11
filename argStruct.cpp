#include "argStruct.h"
#include "shared.h"

#include "io.h"

uint8_t PROGRAM_VERBOSITY_LEVEL = 0;
argStruct* args = NULL;

// TODO:
// check multiple of same argument
// check if given files exist here
argStruct* argStruct_get(int argc, char** argv) {

    if (argc == 0) {
        print_help(stdout);
        IO::outFilesStruct_destroy(outFiles);
        exit(0);
    }

    args = new argStruct;

    while (*argv) {

        char* arv = *argv;
        char* val = *(++argv);

        if (val == NULL) {
            print_help(stdout);
            exit(0);
        }

        // ###########################
        // #    ACTIONS [-doXXXX]    #
        // ###########################
        //
        // Use action commands to specify the action (i.e. analysis) to be performed.
        // Action commands are of the form `-doXXXX <int>`
        //   where XXXX defines the general type of analysis to be performed
        //   and <int> specifies the exact analysis to be performed.
        // For example, `-doAMOVA <int>` specifies to perform AMOVA analysis.
        //   and <int> specifies the exact AMOVA analysis to be performed.
        //   e.g. `-doAMOVA 1` specifies to perform AMOVA analysis with the genotype likelihoods.
        //
        // The following action commands are available:
        //   -doJGTM <int>  : get joint genotypes matrix for each individual pair
        //                    ret: jgtmat
        //   -doDist <int>  : estimate pairwise distance matrix
        //                    req: jgtmat
        //                    ret: dmat
        //   -doAMOVA <int> : perform AMOVA analysis
        //                    req: dmat, metadata, formula
        //   -doEM <int>    : perform EM optimization
        //                    req: DATASOURCE_GL
        //                    ret: jgtmat
        // TODO add doEM 2 for 10gls optim
        //   -doDxy <int>   : estimate Dxy
        //                    req: dmat
        //   -doPhylo <int> : do neighbor-joining
        //                    
        //   -doIbd <int> 	: detect IBD segments
        //
        //   -doMajorMinor <int> : get major and minor alleles for each site

//TODO
        // e.g. `-doAMOVA`  in lut1 -> "args->doAMOVA" in lut2 ->integer in lut3
        // so if `doAMOVA` is used, use multilayer luts to determine the type of required val etc
        if (strcasecmp("-doUnitTests", arv) == 0) {
            args->doUnitTests = atoi(val);
            return(args);
        } else if (strcasecmp("-doAMOVA", arv) == 0) {
            args->doAMOVA = atoi(val);
        } else if (strcasecmp("-doJGTM", arv) == 0) {
            args->doJGTM = atoi(val);
        } else if (strcasecmp("-doEM", arv) == 0) {
            args->doEM = atoi(val);
        } else if (strcasecmp("-doDxy", arv) == 0) {
            args->doDxy = atoi(val);
        } else if ((strcasecmp("-doNeighborJoining", arv) == 0) || (strcasecmp("-doPhylo", arv) == 0)) {
            args->doPhylo = atoi(val);
        } else if ((strcasecmp("--handle-negative-branch", arv) == 0)) {
            args->handle_neg_branch_length = atoi(val);

        } else if (strcasecmp("-doDist", arv) == 0) {
            args->doDist = atoi(val);
        } else if (strcasecmp("-doIbd", arv) == 0) {
            args->doIbd = atoi(val);
        } else if (strcasecmp("-doMajorMinor", arv) == 0) {
            args->doMajorMinor = atoi(val);
        } else if (strcasecmp("--bcf-src", arv) == 0) {
            args->bcfSrc = atoi(val);
        }

        // ################################
        // #    INPUT FILES [--in-XXX]    #
        // ################################
        //
        // Use input file commands to specify the input file types and their filenames.
        // Input file is defined as the file containing the data to be used in the analyses.
        // Input file commands are of the form `--in-XXX <filename>`
        //   where XXX defines the type of input file
        //   and <filename> specifies the filename of the input file.
        // For example, `--in-vcf <filename>` specifies to read the input file as a VCF file.
        //
        // The following input file commands are available:
        //   --in-vcf <filename> : VCF file input
        //   --in-dm <filename> : distance matrix input
        //   --in-dxy <filename> : Dxy file input

        else if ((strcasecmp("--in-vcf", arv) == 0) || (strcasecmp("--input", arv) == 0) || (strcasecmp("-i", arv) == 0)) {
            args->in_vcf_fn = strdup(val);
        } else if ((strcasecmp("--in-dm", arv) == 0)) {
            args->in_dm_fn = strdup(val);
        } else if ((strcasecmp("--in-dxy", arv) == 0)) {
            args->in_dxy_fn = strdup(val);
        } else if ((strcasecmp("--in-majorminor", arv) == 0)) {
            args->in_majorminor_fn = strdup(val);
        } else if ((strcasecmp("--in-ancder", arv) == 0)) {
            args->in_ancder_fn = strdup(val);
        }

        // TODO if -out has full output name and not prefix, detect it
        // e.g. if you were going to write to a file called "out.amova.csv", and "-out out.amova.csv" is given
        //  then extract prefix from that and use as prefix for other outputs as well
        //  this is a useful feature for snakemake etc
        else if ((strcasecmp("--output", arv) == 0) || (strcasecmp("-out", arv) == 0) || (strcasecmp("-o", arv) == 0)) {
            args->out_fnp = strdup(val);
        }

        // #################################################################
        // #    PRINTING COMMANDS                                          #
        // #    [--printXxXxx/--printXX/-pXX <int>]                        #
        // #################################################################
        //
        // Use printing commands to specify the output files to be generated.
        // This is only needed for output files that are not the default output files
        //   associated with the analyses specified by the action commands.
        //
        // Printing commands are of the form `--printXxXxx <int>` or `--printXX <int>` or `-pXX <int>`
        //   where XxXxx is the long form of the file type to be printed,
        //   XX is the short form of the file type (typically  the first letter(s) of the long form),
        //   and <int> defines the compression level of the output file.
        // Output file names are automatically generated based on the value of the `--output` argument
        //   and the file type to be printed.
        //
        // The following printing commands are available:
        //   --printJointGenotypeCountMatrix/-pJGCD <int> : print joint genotype count distribution
        //   --printAmovaTable/-pAT <int> : print AMOVA table
        //   --printDistanceMatrix/-pDM <int> : print distance matrix
        //   --printBlocksTab <0|1> : print tab-delimited blocks file defining the start and end
        //      positions of each block (default: 0 = do not print, 1 = print)
        //
        // The following compression levels are available:
        //   0 : no compression
        //   1 : gzip compression
        //   2 : bgzip compression
        //
        //

        // TODO maybe use hypen style here --print-joint-geno-count-dist to be consistent with other commands
        // TODO maybe make <int> optional and choose a default one, most people won't care about this

        else if ((strcasecmp("--printJointGenotypeCountMatrix", arv) == 0) || (strcasecmp("--printJGCD", arv) == 0) || (strcasecmp("-pJGCD", arv) == 0)) {
            args->printJointGenotypeCountMatrix = atoi(val);
        } else if ((strcasecmp("--printAmovaTable", arv) == 0) || (strcasecmp("--printAT", arv) == 0) || (strcasecmp("-pAT", arv) == 0)) {
            args->printAmovaTable = atoi(val);
        } else if (strcasecmp("--printDistanceMatrix", arv) == 0 || strcasecmp("-pDM", arv) == 0) {
            args->printDistanceMatrix = atoi(val);
        } else if (strcasecmp("--printBlocksTab", arv) == 0) {
            args->printBlocksTab = atoi(val);
        }

        // #######################################################################
        // #    REGION SPECIFICATION COMMANDS                                    #
        // #    [--region/-r] [--regions-tab/-rf/-R] [--regions-bed/-rb/-Rb]    #
        // #######################################################################
        //
        // Use region specification commands to specify the regions to be used in the analyses.
        //
        // There are two types of region specifications:
        //    1) Specifying a single region
        //    2) Specifying a list of regions (with a regions file or a BED file)
        //
        // Region specification commands for a single region are of the form `--region <int>` or `-r <int>`
        //   where <int> is a string of the form `chr:start-end` or `chr:pos` or `chr`.
        //
        // Regions files can be in one of the following formats:
        //   1. BED file
        //      Region specification commands for a list of regions in BED file format
        //        are of the form `--regions-bed <filename>` or `-rb <filename>` or `-Rb <filename>`
        //
        //   2. TAB-delimited genome position file
        //      Region specification commands for a list of regions in regions file format
        //        are of the form `--regions-tab <filename>` or `-rf <filename>` or `-R <filename>`
        //        where <filename> is the name of the file containing the list of regions .
        //
        // -> Both types of regions files should be indexed using the `tabix` program.
        // -> Only one type of region specification command can be used at a time.
        //
        // Warning: Unexpected behavior may occur if the regions file is not in the correct format.
        //   -> The regions file should be sorted by chromosome name and start position.
        //   -> The regions file should not contain any duplicate or overlapping regions.
        //   -> The regions file should not contain any regions that are not present in the VCF file.
        //   -> The regions file should be indexed using the `tabix` program.
        //   -> The regions file should have the same number of columns throughout the file.
        //
        //
        // Coordinate systems
        // -----------------
        // e.g. Chromosome with name "chr1" is of length 8 bases and consists of the following sequence:
        //           A C T G A C T G
        // 0-based   0 1 2 3 4 5 6 7
        // 1-based   1 2 3 4 5 6 7 8
        //
        // - **BED file**
        //     - 0-based
        //     - [start:included, end:excluded)
        //     - Requirements:
        //       - Should be sorted by chromosome name and start position.
        //       - Should be indexed using the `tabix` program.
        //
        // - **TAB-delimited genome position file**
        //     - 1-based
        //     - [start:included, end:included]
        //     - Requirements:
        //       - Should be sorted by chromosome name and start position.
        //       - Should be indexed using the `tabix` program.
        //       - Should have 1, 2, or 3 columns:
        //         - 1 column: <CHR> (chromosome name)
        //		   - 2 columns: <CHR> <POS> (chromosome name and position)
        //		   - 3 columns: <CHR> <START> <END> (chromosome name and start and end positions)
        //
        //     e.g.
        //     "chr1 2"
        //       The second base of the chromosome "chr1" (C)
        //     "chr1 2 3"
        //        The second and third bases of the chromosome "chr1" (C and T)
        //
        // - **--region/-r**
        //     - 1-based
        //     - [start:included, end:included]
        //     - Requirements:
        //       - Should have 1, 2, or 3 columns:
        //         - 1 column: <CHR> (chromosome name)
        //		   - 2 columns: <CHR>:<POS> (chromosome name and position)
        //		   - 3 columns: <CHR>:<START>-<END> (chromosome name and start and end positions)
        //
        //     e.g.
        //     `--region chr1:0-8` == the entire chromosome "chr1" (8 bases)
        //
        //     //TODO check if this coordinate system is the same in angsd regions -r -rf

        else if (strcasecmp("--region", arv) == 0 || strcmp("-r", arv) == 0) {
            args->in_region = strdup(val);
        }

        // N.B. There are two short forms for the regions file command since
        //      -R   easy to remember from bcftools
        //      -rf  people are used to using it in angsd
        else if (strcasecmp("--regions-tab", arv) == 0 || strcasecmp("-rf", arv) == 0 || strcmp("-R", arv) == 0) {
            args->in_regions_tab_fn = strdup(val);
        }

        else if (strcasecmp("--regions-bed", arv) == 0 || strcasecmp("-rb", arv) == 0) {
            args->in_regions_bed_fn = strdup(val);
        }

        // ###################################################################
        // #    BLOCK BOOTSTRAPPING COMMANDS                                 #
        // #    [--block-size/-bs] [--blocks-tab] [--blocks-bed]             #
        // ###################################################################
        //
        // Use block bootstrapping commands to specify the blocks to be used in the block bootstrapping analyses.
        //
        // There are two types of block bootstrapping specifications:
        // - Block size specification
        // - Block list specification
        //
        // Block size specification commands are of the form `--block-size <int>` or `-bs <int>`
        //  where <int> is the size of the blocks. Using the VCF file as input, the blocks are enumerated
        //  by reading the contig sizes from the VCF header and dividing the contigs into blocks of size <int>.
        //
        // Block list specification commands are of the form `--blocks-tab <filename>` or `-bf <filename>`
        //  where <filename> is the name of the file containing the list of blocks.
        //
        // The blocks file should be in one of the following formats:
        //
        //   1. BED file
        //      Region specification commands for a list of regions in BED file format
        //        are of the form `--regions-bed <filename>` or `-rb <filename>` or `-Rb <filename>`
        //
        //
        //   2. TAB-delimited genome position file
        //     - 1-based
        //     - [start:included, end:included]
        //     - Requirements:
        //       - Should be sorted by chromosome name and start position.
        //       - Should be indexed using the `tabix` program.
        //       - Should have 3 columns:
        //		   - <CHR> <START> <END> (chromosome name and start and end positions)

        // -> Requires a VCF file as input data.
        // -> Only one block bootstrapping specification command can be used at a time.

        else if (strcasecmp("--block-size", arv) == 0 || strcasecmp("-bs", arv) == 0) {
            // read block size as float and convert to int
            // to allow for the use of scientific notation (e.g. 1e6)
            args->blockSize = (int)atof(val);
        }

        else if (strcasecmp("--blocks-tab", arv) == 0) {
            args->in_blocks_tab_fn = strdup(val);
        } else if (strcasecmp("--blocks-bed", arv) == 0) {
            args->in_blocks_bed_fn = strdup(val);

            // #####################################
            // #    ARGUMENTS [--long-form/-sf]    #
            // #####################################
            //
            // Use argument commands to specify the parameters to be used in the analyses.
            // Argument commands are of the form `--long-form <int>` or `-sf <int>`
            //   where long form of the argument starts with double dash `--` and separated by hyphen `-`
            //   and short form of the argument starts with single dash `-` and is typically the first letter(s) of the long form
            //

        } else if (strcasecmp("--dxy-groups", arv) == 0) {
            args->dxyGroups = strdup(val);
        } else if (strcasecmp("--dxy-levels", arv) == 0) {
            args->dxyLevels = strdup(val);
        } else if (strcasecmp("--seed", arv) == 0) {
            args->seed = atoi(val);
        } else if (strcasecmp("--metadata", arv) == 0 || strcasecmp("-m", arv) == 0) {
            args->in_mtd_fn = strdup(val);
        } else if ((strcasecmp("--verbose", arv) == 0) || (strcasecmp("-v", arv) == 0)) {
            // allowed range: [0,3]
            int vval = atoi(val);
            if (vval > MAX_PROGRAM_VERBOSITY_LEVEL) {
                vval = MAX_PROGRAM_VERBOSITY_LEVEL;
            }
            PROGRAM_VERBOSITY_LEVEL = vval;
        }

        else if (strcasecmp("--formula", arv) == 0 || strcasecmp("-f", arv) == 0) {
            args->formula = strdup(val);
        }

        // [-sites <filename>]
        //     Filename of the file containing the major and minor alleles for each site in the VCF file.
        //
        //     The file should be in the following format:
        //	   - 1-based indexing for positions [start:included, end:included]
        //	   - Should be sorted by chromosome name and start position.
        //     - Should be tab-delimited.
        //
        // e.g. Chromosome name <TAB> position <TAB> major allele <TAB> minor allele
        //    chr1	100	A	T
        //    chr1	200	A	G
        //
        //    This file can be obtained using ANGSD maf files.


        else if (strcasecmp("-dev", arv) == 0)
            args->printDev = atoi(val);

        else if (strcasecmp("--minInd", arv) == 0)
            args->minInd = atoi(val);
        else if (strcasecmp("--min-af", arv) == 0) {
            args->min_af = atof(val);
        }


        else if (strcasecmp("--ibdseq-errorprop", arv) == 0) {
            args->ibdseq_errorprop = atof(val);
        }

        else if (strcasecmp("--ibdseq-errormax", arv) == 0) {
            args->ibdseq_errormax = atof(val);
        }

        else if (strcasecmp("--ibdseq-minalleles", arv) == 0) {
            args->ibdseq_minalleles = atoi(val);
            if (args->ibdseq_minalleles < 2) {
                ERROR("--ibdseq-minalleles is set to %d. Allowed range of values is 2<value<N", args->ibdseq_minalleles);
            }
        } else if (strcasecmp("--ibdseq-ibdlod", arv) == 0) {
            args->ibdseq_ibdlod = atof(val);
            if (args->ibdseq_ibdlod < 0) {
                ERROR("--ibdseq-ibdlod must be a positive number.");
            }
        } else if (strcasecmp("--ibdseq-ibdtrim", arv) == 0) {
            args->ibdseq_ibdtrim = atof(val);
            if (args->ibdseq_ibdtrim < 0) {
                ERROR("--ibdseq-ibdtrim must be a nonnegative number.");
            }
        }

        else if (strcasecmp("--em-tole", arv) == 0)
            args->tole = atof(val);

        else if ((strcasecmp("-ws", arv) == 0) || (strcasecmp("--windowSize", arv) == 0)) {
            args->windowSize = (int)atof(val);
        }

        else if ((strcasecmp("-nb", arv) == 0) || (strcasecmp("--nBootstraps", arv) == 0)) {
            args->nBootstraps = (int)atof(val);
        }

        else if (strcasecmp("--maxEmIter", arv) == 0) {
            args->maxEmIter = atoi(val);

        } else if (strcasecmp("--nThreads", arv) == 0 || strcasecmp("-P", arv) == 0) {
            args->nThreads = atoi(val);

        } else if (strcasecmp("--rm-invar-sites", arv) == 0) {
            args->rmInvarSites = atoi(val);

        } else if (strcasecmp("-h", arv) == 0 || strcasecmp("--help", arv) == 0) {
            print_help(stdout);
            exit(0);
        } else {
            ERROR("Unknown argument: \'%s\'\n", arv);
        }

        ++argv;
    }

    if (args->doUnitTests) {
        return(args);
    }

    if (NULL == args->out_fnp) {
        args->out_fnp = strdup("amovaput");
        fprintf(stderr, "\n\t-> -out <output_prefix> not set; will use %s as a prefix for output files.\n", args->out_fnp);
    } else {
        fprintf(stderr, "\n\t-> -out <output_prefix> is set to %s. Output files will have this prefix.\n", args->out_fnp);
    }
    IO::outFilesStruct_init(outFiles);


    CHECK_ARG_INTERVAL_INT(args->nThreads, 0, 100, "-P/--nThreads");
    if (0 == args->nThreads) {
        args->nThreads = 1;
    }

    CHECK_ARG_INTERVAL_INT(args->nBootstraps, 0, 100000, "-nb/--nBootstraps");


    args->check_arg_dependencies();

    return args;
}

void argStruct::check_arg_dependencies() {

    if (minInd == 0) {
        LOG("minInd is set to 0; will use sites with data for all individuals.");

    } else if (minInd == -1) {
        LOG("minInd is not set. Default is setting minInd to 2; will use sites that is nonmissing for both individuals in each individual pair.");
        minInd = 2;
    } else if (minInd == 2) {
        LOG("minInd is set to 2; will use sites that is nonmissing for both individuals in a pair.");
    } else if (minInd == 1 || minInd < -1) {
        ERROR("minInd is set to %d. Minimum value allowed for minInd is 2.", minInd);
    }

    if (in_vcf_fn == NULL && in_dm_fn == NULL) {
        fprintf(stderr, "\n[ERROR] Must supply either --in-vcf <VCF_file> or --in-dm <Distance_matrix_file>.\n");
        exit(1);
    } else if (in_vcf_fn != NULL && in_dm_fn != NULL) {
        fprintf(stderr, "\n[ERROR] Cannot use --in-vcf %s and --in-dm %s at the same time.\n", in_vcf_fn, in_dm_fn);
        exit(1);
    }

    if (printDistanceMatrix != 0) {
        fprintf(stderr, "\n[INFO]\t-> -printMatrix %d; will print distance matrix\n", printDistanceMatrix);
    }

    if (seed == -1) {
        seed = time(NULL);
        srand48(seed);
        LOG("Seed is not defined, will use current time as random seed for the random number generator. Seed is now set to: %d.\n", seed);
        WARN("Used the current time as random seed for random number generator. For parallel runs this may cause seed collisions. Hence, it is recommended to set the seed manually using `--seed <INTEGER>`.");
    } else {
        srand48(seed);
        LOG("Seed is set to: %d.\n", seed);
    }

    if (in_dm_fn != NULL) {
        if (doEM != 0) {
            fprintf(stderr, "\n[ERROR]\t-doEM %i cannot be used with -in_dm %s.", doEM, in_dm_fn);
            exit(1);
        }
    }

    if (windowSize != 0) {
        fprintf(stderr, "\n[INFO]\t-> -windowSize %d; will use sliding windows of size %d\n", windowSize, windowSize);
        NEVER;
    }

    //----------------------------------------------------------------------------------//
    // -doDist      defines the method to estimate the pairwise distance matrix
    // (get dmat)
    //              default: 0
    //
    //              0: do not estimate distance matrix
    //TODO
    //              1: Dij
    //              2:
    //              3:


    if (0 == doDist) {
        //
    } else if (1 == doDist) {
        LOG("-doDist %d; will estimate pairwise distance matrix using Dij.", doDist);
    } else {
        ERROR("-doDist %d is not a valid option.", doDist);
    }



    //----------------------------------------------------------------------------------//
    // -doAMOVA
    //             default: 0
    //             0: do not perform AMOVA
    //             1: perform AMOVA 
    //             2: perform AMOVA with block bootstrapping (requires: block size, block file etc)
    //             3: perform AMOVA with permutation test (orig framework, requires: npermut)
    // requires:
    // - formula
    // - metadata
    // - distance matrix

    if (0 == doAMOVA) {
        //
        if (NULL != in_mtd_fn) {
            WARN("-m/--metadata is provided but no analysis requires it; will ignore -m/--metadata.");
        }
        if (NULL != formula) {
            WARN("-f/--formula is provided but no analysis requires it; will ignore -f/--formula.");
        }

    } else if (1 == doAMOVA) {
        IO::requireArgStr(formula, "--formula/-f", "-doAMOVA 1");
        IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doAMOVA 1");

        if (0 == doDist && NULL == in_dm_fn) {
            ERROR("-doAMOVA %d requires either (1) a method to calculate the pairwise distance matrix (-doDist <int>) or (2) a distance matrix file (--in-dm <file>).", doAMOVA);
        }

    } else {
        ERROR("-doAMOVA %d is not a valid option.", doAMOVA);
    }

    if (0 == doEM) {
        //
        if (-1 != tole) {
            ERROR("-tole %e requires -doEM != 0.", tole);
        }
        if (-1 != maxEmIter) {
            ERROR("--maxEmIter %d requires -doEM != 0.", maxEmIter);
        }
    } else if (1 == doEM) {
        if (-1 == tole) {
            tole = 1e-6;
            LOG("-tole is not set, setting to default value %e. Will terminate the EM algorithm if the change in the log-likelihood is less than %e.", tole, tole);
        } else {
            LOG("-tole is set to %e. Will terminate the EM algorithm if the change in the log-likelihood is less than %e.", tole, tole);
        }

        if (-1 == maxEmIter) {
            maxEmIter = 500;
            LOG("-maxEmIter is not set, setting to default value %d. Will terminate the EM algorithm if the number of iterations exceed %d.", maxEmIter, maxEmIter);
        } else {
            LOG("-maxEmIter is set to %d. Will terminate the EM algorithm if the number of iterations exceed %d.", maxEmIter, maxEmIter);
        }

    } else {
        ERROR("-doEM %d is not a valid option.", doEM);
    }

    //----------------------------------------------------------------------------------//
    // -doDxy
    //              default: 0
    //              0: do not estimate dxy
    //              1: estimate dxy from the distance matrix for all groups in all hierarchical levels defined in the metadata file (requires: method to obtain distance matrix, metadata file)

    // : (doDist 1,2 or in_ft == IN_DM), (--metadata-file <METADATA_FILE> or doDxyStr)

    if (0 == doDxy) {
        //
    } else if (1 == doDxy) {
        IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doDxy 1");
        if (args->in_mtd_fn == NULL) {
            ERROR("-doDxy %d requires --metadata/-m <file>.", doDxy);
        }
    } else {
        ERROR("-doDxy %d is not a valid option.", doDxy);
    }

    //----------------------------------------------------------------------------------//
    // -doPhylo
    //              default: 0
    //              0: do not run phylogenetic tree construction
    //              1: construct phylogenetic tree using neighbor joining with individuals as leaf nodes
    //              2: construct phylogenetic tree using neighbor joining with groups as leaf nodes (requires: `-doDxy`)

    if (0 == doPhylo) {
        //
    } else if (1 == doPhylo) {
        // if (0 == doDist && NULL == in_dm_fn) {
        //     ERROR("-doPhylo %d requires a distance matrix (either -doDist <int> or --in-dm <file>).", doPhylo);
        // }
    } else if (2 == doPhylo) {
        if (0 == doDxy && NULL == in_dxy_fn) {
            ERROR("-doPhylo %d requires a dxy matrix (either -doDxy <int> or --in-dxy <file>).", doPhylo);
        }
    } else {
        ERROR("-doPhylo %d is not a valid option.", doPhylo);
    }



    if (args->rmInvarSites) {
    } else {
        //TODO req for some?
        // ERROR("No analyses including invar sites are implemented yet"); //TODO
    }

    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GL) {
        LOG("%s INT+ argument value contains %d. Program will use the GL tag from input VCF file.", "--bcf-src", ARG_INTPLUS_BCFSRC_FMT_GL);
    }
    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GT) {
        LOG("%s INT+ argument value contains %d. Program will use the GT tag from input VCF file.", "--bcf-src", ARG_INTPLUS_BCFSRC_FMT_GT);
    }


    // -------------------------------
    // doMajorMinor
    if (args->doMajorMinor == ARG_DOMAJORMINOR_BCF_REFALT1) {
        LOG("%s is set to %d. Program will use the REF allele as the major allele and the ALT1 allele as the minor allele.", "-doMajorMinor", ARG_DOMAJORMINOR_BCF_REFALT1);
    } else if (args->doMajorMinor == ARG_DOMAJORMINOR_INFILE) {
        if (args->in_majorminor_fn == NULL && args->in_ancder_fn == NULL) {
            ERROR("-doMajorMinor %d requires either --in-majorminor <file> or --in-ancder <file>.", args->doMajorMinor);
        }
        if (args->in_majorminor_fn != NULL && args->in_ancder_fn != NULL) {
            ERROR("Program cannot decide which file to use. -doMajorMinor %d requires either --in-majorminor <file> or --in-ancder <file>, not both.", args->doMajorMinor);
        }
        if (args->in_majorminor_fn != NULL) {
            LOG("%s is set to %d, and %s is provided. Program will use the major and minor alleles from the file %s.", "-doMajorMinor", ARG_DOMAJORMINOR_INFILE, "--in-majorminor", args->in_majorminor_fn);
        }
        if (args->in_ancder_fn != NULL) {
            LOG("%s is set to %d, and %s is provided. Program will use the ancestral and derived alleles from the file %s.", "-doMajorMinor", ARG_DOMAJORMINOR_INFILE, "--in-ancder", args->in_ancder_fn);
        }
    }

    // cases where alleles data is not needed:
    // - em ngl=10
    // - input is not vcf
    if ((args->in_vcf_fn != NULL) && (args->doEM != (ARG_DOEM_10GL))) {
        if (args->doMajorMinor == ARG_DOMAJORMINOR_NONE) {
            ERROR("Program requires major/minor alleles for the specified analyses. Please provide a method to obtain major/minor alleles using -doMajorMinor <int>.");
        }
    }


    // -------------------------------



        //  else if (doPhylo == 1 && doDist == 0 && in_dm_fn == NULL) {
        //     fprintf(stderr, "\n[ERROR]\t-> -doPhylo %i requires -doDist 1 or --in-dm <file>.", doPhylo);
        //     exit(1);
        // }
        // if (doPhylo == 2 && doDxy == 0 && in_dxy_fn == NULL) {
        //     fprintf(stderr, "\n[ERROR]\t-> -doPhylo %i requires -doDxy 1 or --in-dxy <file>.", doPhylo);
        //     exit(1);
        // }

        // TODO using regions with block definitions
        // TODO handle empty blocks
    if (blockSize > 0 && in_blocks_tab_fn != NULL) {
        fprintf(stderr, "\n[ERROR]\t-> `--block-size` cannot be used with `--in-blocks-tab`.\n");
        exit(1);
    }
    if (blockSize > 0 && in_blocks_bed_fn != NULL) {
        fprintf(stderr, "\n[ERROR]\t-> `--block-size` cannot be used with `--in-blocks-bed`.\n");
        exit(1);
    }

    if (in_blocks_tab_fn != NULL || in_blocks_bed_fn != NULL) {
        if (in_regions_tab_fn != NULL || in_regions_bed_fn != NULL || in_region != NULL) {
            fprintf(stderr, "\n[ERROR]\tBlock definitions cannot be used with region definitions, yet.\n");
            exit(1);
        }
    }

    if (args->doEM) {
        if (0 == args->doJGTM) {
            ERROR("-doEM requires -doJGTM");

        }
    }
    if (args->doJGTM) {
        if (0 == args->doEM) {
            // vcf_tags_to_read += 

        }
    }

    // TODO block bootstrapping cannot be used with gt data (bcfSrc==ARG_INTPLUS_BCFSRC_FMT_GT)


    if (PROGRAM_WILL_USE_BCF_FMT_GL && PROGRAM_WILL_USE_BCF_FMT_GT) {
        ERROR("Program cannot use both GL and GT data at the same time.");
    }


}

void argStruct_destroy(argStruct* args) {

    if (args->in_vcf_fn != NULL) {
        FREE(args->in_vcf_fn);
    }

    if (args->in_dm_fn != NULL) {
        FREE(args->in_dm_fn);
    }

    if (args->in_mtd_fn != NULL) {
        FREE(args->in_mtd_fn);
    }

    if (args->in_dxy_fn != NULL) {
        FREE(args->in_dxy_fn);
    }

    if (args->in_region != NULL) {
        FREE(args->in_region);
    }
    if (args->in_regions_tab_fn != NULL) {
        FREE(args->in_regions_tab_fn);
    }

    if (args->in_regions_bed_fn != NULL) {
        FREE(args->in_regions_bed_fn);
    }

    if (args->in_blocks_tab_fn != NULL) {
        FREE(args->in_blocks_tab_fn);
    }

    if (args->in_blocks_bed_fn != NULL) {
        FREE(args->in_blocks_bed_fn);
    }

    if (args->in_majorminor_fn != NULL) {
        FREE(args->in_majorminor_fn);
    }

    if (args->in_ancder_fn != NULL) {
        FREE(args->in_ancder_fn);
    }

    FREE(args->out_fnp);

    if (args->formula != NULL) {
        FREE(args->formula);
    }

    if (args->keyCols != NULL) {
        FREE(args->keyCols);
    }
    if (args->doDxyStr != NULL) {
        FREE(args->doDxyStr);
    }

    delete args;
}