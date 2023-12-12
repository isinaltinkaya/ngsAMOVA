#include "argStruct.h"

#include "io.h"

//  default: 0 (verbose mode off)
u_char VERBOSE = 0;

argStruct *args = NULL;

// TODO check multiple of same argument
/// @brief argStruct_get read command line arguments
/// @param argc
/// @param argv
/// @return pointer to argStruct
argStruct *argStruct_get(int argc, char **argv) {
    args = new argStruct;

    while (*argv) {
        // TODO check if given files exist here

        char *arv = *argv;
        char *val = *(++argv);

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
        //   -doAMOVA <int> : perform AMOVA analysis
        //   -doEM <int> : perform EM optimization
        //   -doDxy <int> : estimate Dxy
        //   -doPhylo <int> : do neighbor-joining
        //   -doDist <int> : estimate pairwise distance matrix
		//   -doIbd <int> 	: detect IBD segments

        // TODO idea multilayer argument reading using LUTs
        // e.g. `-doAMOVA`  in lut1 -> "args->doAMOVA" in lut2 ->integer in lut3
        // so if `doAMOVA` is used, use multilayer luts to determine the type of required val etc
        // getArg("-doAMOVA", arv);

        if (strcasecmp("-doAMOVA", arv) == 0) {
            args->doAMOVA = atoi(val);
        } else if (strcasecmp("-doEM", arv) == 0) {
            args->doEM = atoi(val);
        } else if (strcasecmp("-doDxy", arv) == 0) {
            if (strIsNumeric(val)) {
                args->doDxy = atoi(val);
            } else {
                args->doDxy = 999;  // indicates that a string is provided
                args->doDxyStr = strdup(val);
            }
        } else if ((strcasecmp("-doNeighborJoining", arv) == 0) || (strcasecmp("-doPhylo", arv) == 0)) {
            args->doPhylo = atoi(val);
        } else if ((strcasecmp("--handle-negative-branch", arv) == 0) ){
		args->handle_neg_branch_length=atoi(val);

        } else if (strcasecmp("-doDist", arv) == 0) {
            args->doDist = atoi(val);
        } else if (strcasecmp("-doIbd", arv) == 0) {
            args->doIbd = atoi(val);
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
            if (atoi(val) == 0) {
                // explicit verbose off, use the default value 0
                fprintf(stderr, "[INFO]\t-> Verbosity disabled explicitly. Will not print any information.\n");
            } else if (strIsNumeric(val)) {
                int verbose_val = atoi(val) - 1;
                if (verbose_val > MAX_VERBOSE_LEVEL) {
                    verbose_val = MAX_VERBOSE_LEVEL;
                    fprintf(stderr, "\n[WARNING]\tVerbosity level is set to %d, which is greater than the maximum verbosity level %d. Setting verbosity level to MAX.\n", verbose_val, MAX_VERBOSE_LEVEL);
                }
                BITSET(VERBOSE, verbose_val);
            } else {
                ERROR("Invalid verbosity level provided: %s\n", val);
            }
            fprintf(stderr, "\n[INFO]\t-> Verbosity level is set to %d.\n", WHICH_BIT_SET1(VERBOSE));
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
        else if (strcasecmp("--min-af", arv) == 0){
            args->min_af = atof(val);
		}


        else if (strcasecmp("--ibdseq-errorprop", arv) == 0){
			args->ibdseq_errorprop = atof(val);
		}

        else if (strcasecmp("--ibdseq-errormax", arv) == 0){
			args->ibdseq_errormax = atof(val);
		}

        else if (strcasecmp("--ibdseq-minalleles", arv) == 0){
            args->ibdseq_minalleles = atoi(val);
			if(args->ibdseq_minalleles<2){
				ERROR("--ibdseq-minalleles is set to %d. Allowed range of values is 2<value<N");
			}
		}
        else if (strcasecmp("--ibdseq-ibdlod", arv) == 0){
            args->ibdseq_ibdlod= atof(val);
			if(args->ibdseq_ibdlod<0){
				ERROR("--ibdseq-ibdlod must be a positive number.");
			}
		}
        else if (strcasecmp("--ibdseq-ibdtrim", arv) == 0){
            args->ibdseq_ibdtrim= atof(val);
			if(args->ibdseq_ibdtrim<0){
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

        else if (strcasecmp("--maxEmIter", arv) == 0)
            args->maxEmIter = atoi(val);

        else if (strcasecmp("--isSim", arv) == 0) {
            args->isSim = atoi(val);

        } else if (strcasecmp("--nThreads", arv) == 0 || strcasecmp("-P", arv) == 0) {
            args->nThreads = atoi(val);

        } else if (strcasecmp("-h", arv) == 0 || strcasecmp("--help", arv) == 0) {
            print_help(stdout);
            exit(0);
        } else {
            ERROR("Unknown argument: \'%s\'\n", arv);
        }

        ++argv;
    }

    if (NULL == args->out_fnp) {
        args->out_fnp = strdup("amovaput");
        fprintf(stderr, "\n\t-> -out <output_prefix> not set; will use %s as a prefix for output files.\n", args->out_fnp);
    } else {
        fprintf(stderr, "\n\t-> -out <output_prefix> is set to %s. Output files will have this prefix.\n", args->out_fnp);
    }
    IO::outFilesStruct_init(outFiles);

    args->check_arg_dependencies();
    args->print();
    args->print(stderr);

    return args;
}

void argStruct::check_arg_dependencies() {
    if (0 != printJointGenotypeCountMatrix) {
        // TODO implement this
        ERROR("printJointGenotypeCountMatrix is not available.");
    }

    if (minInd == 0) {
        ARGLOG("minInd is set to 0; will use sites with data for all individuals.")

    } else if (minInd == -1) {
        ARGLOG("minInd is not set. Default is setting minInd to 2; will use sites that is nonmissing for both individuals in each individual pair.")
        minInd = 2;
    } else if (minInd == 2) {
        ARGLOG("minInd is set to 2; will use sites that is nonmissing for both individuals in a pair.")
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

    // if data is from in_dm_fn (distance matrix file) can only use -doDist 3
    if (NULL != in_dm_fn && doDist != 3) {
        ERROR("--in-dm <DM_file> requires -doDist 3.");
    }

    if (seed == -1) {
        seed = time(NULL);
        srand48(seed);
        ARGLOG("Seed is not defined, will use current time as random seed for the random number generator. Seed is now set to: %d.\n", seed);
        WARNING("Used the current time as random seed for random number generator. For parallel runs this may cause seed collisions. Hence, it is recommended to set the seed manually using `--seed <INTEGER>`.");
    } else {
        srand48(seed);
        ARGLOG("Seed is set to: %d.\n", seed);
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
    // -doDist      defines the method to estimate distance matrix
    //              default: 0
    //
    //              0: do not estimate distance matrix
    //              1: estimate distance matrix from genotype likelihoods
    //              2: estimate distance matrix from genotypes
    //              3: read distance matrix from file

	//TODO
    // if (0 == doDist) {
        // fprintf(stderr, "\n\n-doDist 0: Nothing to do; will exit.\n\n");
        // exit(0);
    // } else
		if (1 == doDist) {
        LOG("-doDist 1, will estimate distance matrix from genotype likelihoods.");

        IO::requireArgFile(in_vcf_fn, "--in-vcf", "-doDist 1");

        if (0 == doEM) {
            ERROR("-doDist 1 requires -doEM != 0.\n");
        }
        squareDistance = 1;
    } else if (2 == doDist) {
        IO::requireArgFile(in_vcf_fn, "--in-vcf", "-doDist 2");
        squareDistance = 1;
    } else if (3 == doDist) {
        IO::requireArgFile(in_dm_fn, "--in-dm", "-doDist 3");
        squareDistance = 0;

        if (0 != printDistanceMatrix) {
            ERROR("-doDist 3 cannot be used with --printDistanceMatrix.\n");
        }
	}
		//TODO
    // } else {
        // ERROR("-doDist %d is not a valid option.", doDist);
    // }

    //----------------------------------------------------------------------------------//
    // -doAMOVA
    //             default: 0
    //             0: do not perform AMOVA
    //             1: perform AMOVA

    if (0 == doAMOVA) {
        //
        if (NULL != in_mtd_fn) {
            WARNING("-m/--metadata is provided but no analysis requires it; will ignore -m/--metadata.");
        }

        if (NULL != formula) {
            WARNING("-f/--formula is provided but no analysis requires it; will ignore -f/--formula.");
        }

    } else if (1 == doAMOVA) {
        IO::requireArgStr(formula, "--formula/-f", "-doAMOVA 1");
        IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doAMOVA 1");

        if (0 == doDist) {
            ERROR("-doAMOVA %d requires -doDist != 0.", doAMOVA);
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
    //              1: estimate dxy from the distance matrix for all groups in each hierarchical level defined in the metadata file (requires: method to obtain distance matrix, metadata file)
    //              2: estimate dxy from the distance matrix for the groups in metadata specified by the user via --dxy-groups (requires: method to obtain distance matrix, metadata file, --dxy-groups)
    //              3: estimate dxy from the distance matrix for the groups in the hierarchical level specified by the user via --dxy-levels (requires: method to obtain distance matrix, metadata file, --dxy-levels)

    // : (doDist 1,2 or in_ft == IN_DM), (--metadata-file <METADATA_FILE> or doDxyStr)
    if (0 == doDxy) {
        //
    } else if (1 == doDxy) {
        IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doDxy 1");
        if (0 == doDist) {
            ERROR("-doDxy %d requires -doDist != 0.", doDxy);
        }
    } else if (2 == doDxy) {
        IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doDxy 2");
        IO::requireArgStr(dxyGroups, "--dxy-groups", "-doDxy 2");
        if (0 == doDist) {
            ERROR("-doDxy %d requires -doDist != 0.", doDxy);
        }
    } else if (3 == doDxy) {
        IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doDxy 3");
        IO::requireArgStr(dxyLevels, "--dxy-levels", "-doDxy 3");
        if (0 == doDist) {
            ERROR("-doDxy %d requires -doDist != 0.", doDxy);
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

    if (0 == doPhylo){
        //
    } else if (1 == doPhylo) {
        if (0 == doDist && NULL == in_dm_fn) {
            ERROR("-doPhylo %d requires a distance matrix (either -doDist <int> or --in-dm <file>).", doPhylo);
        }
    } else if (2 == doPhylo) {
        if (0 == doDxy && NULL == in_dxy_fn) {
            ERROR("-doPhylo %d requires a dxy matrix (either -doDxy <int> or --in-dxy <file>).", doPhylo);
        }
    } else {
        ERROR("-doPhylo %d is not a valid option.", doPhylo);
    }

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
}

void argStruct_destroy(argStruct *args) {
    FFREE(args->in_vcf_fn);
    FFREE(args->in_dm_fn);
    FFREE(args->in_mtd_fn);
    FFREE(args->in_dxy_fn);

    FFREE(args->in_region);
    FFREE(args->in_regions_tab_fn);
    FFREE(args->in_regions_bed_fn);
    FFREE(args->in_blocks_tab_fn);
    FFREE(args->in_blocks_bed_fn);

    FFREE(args->out_fnp);

    FFREE(args->formula);
    FFREE(args->keyCols);
    FFREE(args->doDxyStr);

    // check below

    delete args;
}

// TODO
void argStruct::print() {
    if (NULL == outFiles->out_args_fs->kbuf) {
        outFiles->out_args_fs->kbuf = kbuf_init();
    }
    kstring_t *kbuf = outFiles->out_args_fs->kbuf;

    kputs("\nCommand: ", kbuf);
    // ksprintf(kbuf, "%s", args->command);
}

void argStruct::print(FILE *fp) {
    // fprintf(fp, "\nCommand: %s", args->command);//TODO
    fprintf(fp, "\n");
}
