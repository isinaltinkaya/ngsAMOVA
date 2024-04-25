#include "argStruct.h"
#include "shared.h"
#include "version.h"
#include "build.h"

#include <htslib/hts.h> // hts_version()

#include "io.h"

uint8_t PROGRAM_VERBOSITY_LEVEL = 0;
char* PROGRAM_VERSION_INFO = NULL;
char* PROGRAM_COMMAND = NULL;
argStruct* args = NULL;




void version_page() {
    fprintf(stderr, "ngsAMOVA [version: %s] [build: %s %s] [htslib: %s]\n", NGSAMOVA_VERSION, __DATE__, __TIME__, hts_version());
    fprintf(stderr, "\n");
    fprintf(stderr, "Build details:\n");
    fprintf(stderr, "\t-> CXX=%s\n", NGSAMOVA_MAKE_CXX);
    fprintf(stderr, "\t-> CXXFLAGS=%s\n", NGSAMOVA_MAKE_CXXFLAGS);
    fprintf(stderr, "\t-> CPPFLAGS=%s\n", NGSAMOVA_MAKE_CPPFLAGS);
    fprintf(stderr, "\t-> LIBS=%s\n", NGSAMOVA_MAKE_LIBS);
    fprintf(stderr, "\t-> FLAGS=%s\n", NGSAMOVA_MAKE_FLAGS);
    fprintf(stderr, "\t-> HTSSRC=%s\n", NGSAMOVA_MAKE_HTSSRC);
    fprintf(stderr, "\n");

}



// TODO:
// check multiple of same argument
// check if given files exist here
argStruct* argStruct_get(int argc, char** argv) {

    if (argc == 0) {
        print_help(stdout);
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
        // For example, `-doAMOVA <int>` specifies to perform AMOVA analysis.
        //
        // The following action commands are available:
        //   -doIbd <int> 	: detect IBD segments
        //


// --in-vcf <STR> : input VCF file
//                 returns: vcfdata, gldata

// --in-dm <STR>  : input distance matrix file
//                 returns: dmat

// -doJGTM <INT>  : get pairwise joint genotype matrix
//                 returns: jgtmat

// -doDist <INT>  : estimate pairwise distance matrix
//                 requires: jgtmat
//                 returns: dmat

// -doAMOVA <INT> : perform AMOVA analysis
//                 requires: dmat, metadata, formula

// -doEM <INT>    : perform EM optimization
//                 requires: gldata
//                 returns: jgtmat

// -doDxy <INT>   : estimate Dxy
//                 requires: dmat

// -doPhylo <INT> : do neighbor-joining tree
//                 requires: dmat

// -doMajorMinor <INT> : get major and minor alleles for each site 
//                 requires: vcfdata, alleles input file (optional)
//                 returns: alleles



        if (strcasecmp("-doUnitTests", arv) == 0) {
            args->doUnitTests = atoi(val);
            return(args);
        } else if (strcasecmp("-doAMOVA", arv) == 0) {
            args->doAMOVA = atoi(val);
        } else if (strcasecmp("-doJGTM", arv) == 0) {
            args->doJGTM = atoi(val);
        } else if (strcasecmp("-doEM", arv) == 0) {
            args->doEM = atoi(val);
            // when performing block bootstrapping em, use --verbose >1 to print per em optimization termination reason, last d value, and number of iterations
            // TODO add this to help
        } else if (strcasecmp("-doDxy", arv) == 0) {
            args->doDxy = atoi(val);
        } else if ((strcasecmp("-doNeighborJoining", arv) == 0) || (strcasecmp("-doPhylo", arv) == 0)) {
            args->doPhylo = atoi(val);
        } else if ((strcasecmp("--handle-negative-branch", arv) == 0)) {
            args->handle_neg_branch_length = atoi(val);

        } else if (strcasecmp("-doDist", arv) == 0) {
            args->doDist = atoi(val);
        } else if (strcasecmp("--dm-method", arv) == 0) {
            args->dm_method = atoi(val);
        } else if (strcasecmp("--dm-transform", arv) == 0) {
            args->dm_transform = atoi(val);
        } else if (strcasecmp("-doIbd", arv) == 0) {
            args->doIbd = atoi(val);
        } else if (strcasecmp("-doMajorMinor", arv) == 0) {
            args->doMajorMinor = atoi(val);
        } else if (strcasecmp("-doBlockBootstrap", arv) == 0) {
            args->doBlockBootstrap = atoi(val);
        } else if (strcasecmp("--bcf-src", arv) == 0) {
            args->bcfSrc = atoi(val);
        } else if (strcasecmp("--prune-dmat", arv) == 0) {
            // any downstream analysis will use the pruned version of the dmat
            // --print-dm 1 will print the original dmat only
            // --print-dm 2 will print the pruned dmat only
            // --print-dm 3 will print both
            // see ARG_INTPLUS_PRINT_DM_*
            args->prune_dmat = atoi(val);
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
        //   --print-jgtm <int> : print joint genotype matrix
        //             requires : --in-vcf, -doJGTM
        //   --print-dm <int> : print distance matrix
        //             requires : -doDist

        //   --print_amova/-pAT <int> : print AMOVA table
        //   --print_blocks <0|1> : print tab-delimited blocks file defining the start and end
        //      positions of each block (default: 0 = do not print, 1 = print)
        //
        // The following compression levels are available:
        //   0 : no compression
        //   1 : gzip compression
        //   2 : bgzip compression
        //
        //

        // TODO maybe use hypen style here --print-joint-geno-count-dist to be consistent with other commands

        else if (strcasecmp("--print-jgtm", arv) == 0) {
            args->print_jgtm = atoi(val);
        } else if (strcasecmp("--print-dm", arv) == 0) {
            args->print_dm = atoi(val);
        } else if ((strcasecmp("--print-amova", arv) == 0)) {
            args->print_amova = atoi(val);
        } else if (strcasecmp("--print-blocks", arv) == 0) {
            args->print_blocks = atoi(val);
        } else if (strcasecmp("--print-bootstrap", arv) == 0) {
            // prints to <prefix>.bootstrap_samples.tsv
            // with header:
            // Rep\tPos\tBlockID\tBlockContig\tBlockStart\tBlockEnd
            // Rep 	Pos 	BlockID 	BlockContig BlockStart 	BlockEnd
            // Rep:        Replicate number
            // Pos:        Position of the sampled block in the replicates synthetic genome (0-based)
            // BlockID:    ID of the sampled block
            // BlockContig: Name of the chromosome to which the block belongs to
            // BlockStart: 1-based, inclusive [start position of the block with the given BlockID
            // BlockEnd:   1-based, inclusive end] position of the block with the given BlockID
            //
            // e.g.
            // if we have 4 blocks in the original genome, each with a size of 1000 bps 
            // {block1, block2, block3, block4} 
            // (0     , 1     , 2     , 3) BlockIDs
            // internal representation start end pos: (0-based inclusive start exclusive end)
            // (all belonging to chr1)
            // block1 start:0,    end:1000 
            // block2 start:1000, end:2000
            // block3 start:2000, end:3000
            // block4 start:3000, end:4000
            // blocks tab representation: (1-based inclusive start end)
            // block1 start:1,    end:1000
            // block2 start:1001, end:2000
            // block3 start:2001, end:3000
            // block4 start:3001, end:4000
            //
            // and we have one replicate which has sampled the following ordered set:
            // {block3, block1, block2, block4}
            //
            // Our output file output.bootstrap_samples.tsv will have:
            // Rep 	Pos 	BlockID 	BlockContig BlockStart 	BlockEnd
            // 0    0       2           chr1        2001        3000
            // 0    1       0           chr1        0           1000
            // 0    2       1           chr1        1001        2000
            // 0    3       3           chr1        3001        4000
            // 

            args->print_bootstrap = atoi(val);
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

        else if ((strcasecmp("--windowSize", arv) == 0)) {
            args->windowSize = (int)atof(val);
        }

        else if ((strcasecmp("--nBootstraps", arv) == 0)) {
            args->nBootstraps = (int)atof(val);
        }

        else if ((strcasecmp("--bootstrap-ci", arv) == 0)) {
            args->bootstrap_pctci = atof(val);
        }

        else if ((strcasecmp("--maxEmIter", arv) == 0) || (strcasecmp("--em-max-iter", arv) == 0)) {
            args->maxEmIter = atoi(val);

        } else if ((strcasecmp("-P", arv) == 0) || (strcasecmp("-@", arv) == 0) || (strcasecmp("--nThreads", arv) == 0) || (strcasecmp("--threads", arv) == 0) || (strcasecmp("-nThreads", arv) == 0)) {
            args->nThreads = atoi(val);

        } else if (strcasecmp("--rm-invar-sites", arv) == 0) {
            args->rmInvarSites = atoi(val);

        } else if (strcasecmp("--drop-pairs", arv) == 0) {
            // if an individual pair has no shared sites, should the program drop the pair or give err and exit?
            args->allow_mispairs = atoi(val);

        } else if (strcasecmp("--min-pairsites", arv) == 0) {
            // TODO
            // drop the pair if the number of shared sites is less than this value
            args->pair_min_n_sites = (int)atof(val);

        } else if (strcasecmp("--min-npairs", arv) == 0) {
            args->min_n_pairs = atoi(val);

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


    if (-1 != args->allow_mispairs) {
        CHECK_ARG_INTERVAL_01(args->allow_mispairs, "--drop-pairs");
        WARN("--drop pairs is set to %d. For downstream analyses, use this option with caution as it may lead to decreased power.", args->allow_mispairs);
    } else {
        args->allow_mispairs = 0;
        LOG("--drop-pairs is not set, setting to default value %d (do not drop pairs). Program will give error and exit if an individual pair has no shared sites.", args->allow_mispairs);
    }

    if (-1 != args->pair_min_n_sites) {
        CHECK_ARG_INTERVAL_INT(args->pair_min_n_sites, 1, INT_MAX, "--min-pairsites");
    } else {
        args->pair_min_n_sites = 1;
        LOG("--min-pairsites is not set, setting to default value %d (perform the action defined via --drop-pairs if an individual pair has 0 shared sites).", args->pair_min_n_sites);
    }

    if (args->print_dm & ARG_INTPLUS_PRINT_DM_PRUNED) {
        if (args->prune_dmat != 1) {
            ERROR("Cannot print pruned distance matrix when pruning is not enabled.");
        }
    }


    if (-1 != args->min_n_pairs) {
        // requires: --drop-pairs 1
        CHECK_ARG_INTERVAL_INT(args->min_n_pairs, 0, INT_MAX, "--min-npairs");
        // if (0 == args->allow_mispairs) {
        //     ERROR("Cannot set --min-npairs without setting --drop-pairs to 1. Please set --drop-pairs to 1.");
        // }
    } else {
        args->min_n_pairs = 1;
        LOG("--min-npairs is not set, setting to default value %d. Program will give error if the number of individual pairs left after dropping pairs is less than this value.", args->min_n_pairs);
    }

    CHECK_ARG_INTERVAL_INT(args->nThreads, 0, 100, "-P/--nThreads");
    if (0 == args->nThreads) {
        args->nThreads = 1;
    }

    if (NULL != args->in_vcf_fn) {
        fprintf(stderr, "\n[INFO]\tFound input VCF file: %s\n", args->in_vcf_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_VCF;
        if (-1 == args->bcfSrc) {
            ERROR("BCF data source is necessary for VCF input. Please set the BCF data source using --bcf-src.");
        }
    }

    if (NULL != args->in_dm_fn) {
        fprintf(stderr, "\n[INFO]\tFound input distance matrix file: %s\n", args->in_dm_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_DM;
    }

    if (NULL != args->in_dxy_fn) {
        //TODO ??
        fprintf(stderr, "\n[INFO]\tFound input dxy file: %s\n", args->in_dxy_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_DXY;
    }

    if (NULL != args->in_mtd_fn) {
        fprintf(stderr, "\n[INFO]\tFound input metadata file: %s\n", args->in_mtd_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_METADATA;
    }

    if (NULL != args->in_majorminor_fn) {
        fprintf(stderr, "\n[INFO]\tFound input major/minor alleles file: %s\n", args->in_majorminor_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_MAJORMINOR;
    }

    if (NULL != args->in_ancder_fn) {
        fprintf(stderr, "\n[INFO]\tFound input ancestral/derived alleles file: %s\n", args->in_ancder_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_ANCDER;
    }

    if (NULL != args->in_blocks_bed_fn) {
        fprintf(stderr, "\n[INFO]\tFound input blocks BED file: %s\n", args->in_blocks_bed_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_BLOCKS;
    }

    if (NULL != args->in_blocks_tab_fn) {
        fprintf(stderr, "\n[INFO]\tFound input blocks TSV file: %s\n", args->in_blocks_tab_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_BLOCKS;
    }

    if (NULL != args->in_regions_bed_fn) {
        fprintf(stderr, "\n[INFO]\tFound input regions BED file: %s\n", args->in_regions_bed_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_REGIONS;
    }

    if (NULL != args->in_regions_tab_fn) {
        //TODO test
        fprintf(stderr, "\n[INFO]\tFound input regions TSV file: %s\n", args->in_regions_tab_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_REGIONS;
    }

    if ((PROGRAM_HAS_INPUT_DM || PROGRAM_HAS_INPUT_MULTIDM) && PROGRAM_HAS_INPUT_VCF) {
        ERROR("Both VCF and distance matrix input are provided. Only one type of input is allowed at the same time.");
    }

    if (PROGRAM_HAS_INPUT_DM) {
        if (args->blockSize > 0) {
            ERROR("-blockSize is not supported for distance matrix input.");
        }
        if (args->doEM) {
            ERROR("-doEM is not available for distance matrix input.");
        }
        if (args->doJGTM) {
            ERROR("-doJGTM is not available for distance matrix input.");
        }
    }

    if (PROGRAM_HAS_INPUT_METADATA && (!PROGRAM_NEEDS_METADATA)) {
        WARN("Metadata file is provided but no analysis requires it; will ignore the metadata file %s.", args->in_mtd_fn);
        if (NULL != args->formula) {
            WARN("Formula is set but no analysis requires it; will ignore the formula '%s'.", args->formula);
        }
    }

    if ((!(PROGRAM_HAS_INPUT_METADATA)) && PROGRAM_NEEDS_METADATA) {
        ERROR("Metadata file is required for the specified analyses. Please provide the metadata file using --metadata.");
    }

    if (PROGRAM_HAS_INPUT_METADATA && PROGRAM_NEEDS_METADATA) {
        if (args->formula == NULL) {
            ERROR("Metadata file is provided but no formula is set. Please set the formula using -f/--formula.");
        }
    }


    if (args->doBlockBootstrap){
        if (PROGRAM_WILL_USE_BCF_FMT_GT) {
            ERROR("Program cannot use GT data for block bootstrapping.");
        }


        if (args->nBootstraps != -1) {
            CHECK_ARG_INTERVAL_INT(args->nBootstraps, 1, 100000, "-nb/--nBootstraps");
        }

        if (-1.0 == args->bootstrap_pctci) {
            args->bootstrap_pctci = 95.0;
            LOG("Bootstrap confidence interval is not set, setting to default value %.17g %%.", args->bootstrap_pctci);

        } else {
            CHECK_ARG_INTERVAL_EE_DBL(args->bootstrap_pctci, 0.0, 100.0, "--bootstrap-ci");
        }

        if (args->blockSize == 0 && args->in_blocks_tab_fn == NULL && args->in_blocks_bed_fn == NULL) {
            ERROR("Block bootstrapping is enabled but no block specification is provided. Please provide block size or block file.");
        }
        if (args->nBootstraps == -1) {
            ERROR("Block bootstrapping is enabled but the number of bootstraps is not set. Please set the number of bootstraps using -nb/--nBootstraps.");
        }
    }

    if (PROGRAM_HAS_INPUT_METADATA && (!PROGRAM_NEEDS_METADATA)) {
        ERROR("Metadata file is provided but none of the specified analyses use it. Please remove the metadata file or set the correct analysis.");

    }


    if (args->doAMOVA) {
        if (ARG_INTPLUS_UNSET == args->print_amova) {
            args->print_amova = ARG_INTPLUS_PRINT_AMOVA_CSV;
            LOG("--print-amova is not set, setting to default value %d (print AMOVA table in CSV format).", args->print_amova);
        }
    }

    args->check_arg_dependencies();


    ASSERT(asprintf(&PROGRAM_VERSION_INFO, "ngsAMOVA [version: %s] [build: %s %s] [htslib: %s]", NGSAMOVA_VERSION, __DATE__, __TIME__, hts_version()) != -1);

    ASSERT(asprintf(&PROGRAM_COMMAND, "%s", argv[0]) != -1);


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
    //              1: calculate the pairwise distance matrix using the method defined via --dm-method
    //              2: read the distance matrix from the file defined via --in-dm
    //              3: do doDist 1 for original dataset and the block bootstrapped datasets

    if (0 == doDist) {
        //
    } else if ((1 == doDist) || (3 == doDist)) {
        // if dm_method dm_transform not defined, set to default and inform user
        if (-1 == dm_method) {
            dm_method = DMAT_METHOD_DIJ;
            LOG("Distance matrix method is not set, setting to default value %d (Dij).", dm_method);
        } else {
            CHECK_ARG_INTERVAL_INT(dm_method, 0, 9, "--dm-method");
            LOG("--dm-method is set to %d (%s).", dm_method, get_dmat_method_str(dm_method));
        }
        if (-1 == dm_transform) {
            dm_transform = DMAT_TRANSFORM_NONE;
            LOG("Distance matrix transform is not set, setting to default value %d (none).", dm_transform);
        } else {
            CHECK_ARG_INTERVAL_INT(dm_transform, 0, 1, "--dm-transform");
        }

        if (ARG_INTPLUS_UNSET == args->print_dm) {
            args->print_dm = ARG_INTPLUS_PRINT_DM_ORIG;
            LOG("--print-dm is not set, setting to default value %d (print the distance matrix).", args->print_dm);
        }

    } else if (2 == doDist) {
        if (NULL == in_dm_fn) {
            ERROR("-doDist %d requires a distance matrix file (--in-dm <file>).", doDist);
        }
        LOG("-doDist %d; will read the distance matrix from the file %s.", doDist, in_dm_fn);
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

    if (doAMOVA) {
        if (ARG_DOAMOVA_SINGLERUN == doAMOVA) {
            IO::requireArgStr(formula, "--formula/-f", "-doAMOVA 1");
            IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doAMOVA 1");

            if (0 == doDist && NULL == in_dm_fn) {
                ERROR("-doAMOVA %d requires either (1) a method to calculate the pairwise distance matrix (-doDist <int>) or (2) a distance matrix file (--in-dm <file>).", doAMOVA);
            }

        } else if (ARG_DOAMOVA_BOOTSTRAP == doAMOVA) {
            IO::requireArgStr(formula, "--formula/-f", "-doAMOVA 1");
            IO::requireArgFile(in_mtd_fn, "--metadata/-m", "-doAMOVA 1");

            if (0 == doDist && NULL == in_dm_fn) {
                ERROR("-doAMOVA %d requires either (1) a method to calculate the pairwise distance matrix (-doDist <int>) or (2) a distance matrix file (--in-dm <file>).", doAMOVA);
            }
            if (!PROGRAM_WILL_USE_BCF_FMT_GL) {
                ERROR("-doAMOVA %d requires genotype likelihoods in the input VCF file (--bcf-src %d).", doAMOVA, ARG_INTPLUS_BCFSRC_FMT_GL);
            }
        } else if (ARG_DOAMOVA_PERMTEST) {
            NEVER;//TODO 
        } else {
            ERROR("-doAMOVA %d is not a valid option.", doAMOVA);
        }
    }

    // If the input data source is genotype likelihoods(`- - bcf - src 1`), the EM optimization(`-doEM 1`) is needed to obtain the joint genotype matrix.
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

    // : (doDist 1,2 or in_ft == IN_DM), (--metadata-file <METADATA_FILE>)

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
    } else {
        ERROR("-doPhylo %d is not a valid option.", doPhylo);
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
        if (args->doMajorMinor == ARG_DOMAJORMINOR_UNSET) {
            ERROR("Program requires major/minor alleles for the specified analyses. Please provide a method to obtain major/minor alleles using -doMajorMinor <int>.");
        }
    }


    // -------------------------------

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

    if (PROGRAM_WILL_USE_BCF_FMT_GL && PROGRAM_WILL_USE_BCF_FMT_GT) {
        ERROR("Program cannot use both GL and GT data at the same time.");
    }

    if (args->rmInvarSites && PROGRAM_WILL_USE_BCF_FMT_GL) {
        LOG("(--rm-invar-sites %d --bcf-src %d) Program will remove invariant sites based on AD tag from input VCF file.", args->rmInvarSites, args->bcfSrc);
    }

    if (args->rmInvarSites && PROGRAM_WILL_USE_BCF_FMT_GL && PROGRAM_WILL_USE_BCF_FMT_GT) {
        ERROR("Program cannot remove invariant sites when using both GL and GT data.");
    }
    // todo use rmInvarSites INT value to determine what to use to remove invar

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

    FREE(PROGRAM_VERSION_INFO);
    FREE(PROGRAM_COMMAND);

    delete args;
}


/// @brief usage - print usage
/// @param fp pointer to the file to print to
void print_help(FILE* fp) {
    fprintf(fp,
        "\n"
        "Program: ngsAMOVA\n");

    fprintf(fp,
        "\n"
        "Usage:\tngsAMOVA <command> [options]\n"
        "\n"
        "Tool for performing the Analysis of Molecular Variance [AMOVA]\n"
        "\n"
        "\n");

    // fprintf(fp, "\n");
    // fprintf(fp, "Usage: ngsAMOVA [options] -i <vcf file> -out <output file>\n");
    // fprintf(fp, "\n");
    // fprintf(fp, "Options:\n");
    // fprintf(fp, "  -in <vcf file>		: input vcf file\n");
    // fprintf(fp, "  -out <output file>		: output file\n");
    // fprintf(fp, "  -doAMOVA <0/1>		: do AMOVA (default: 1)\n");
    // fprintf(fp, "  -doTest <0/1>		: do EM test (default:
    // fprintf(fp, "  -printMatrix <0/1>		: print distance matrix (default: 0)\n");
    // fprintf(fp, "  -printSFS <0/1>		: print SFS (default: 0)\n");

    //
    //
}

/// AMOVA Formula Specification:
/// ----------------------------
/// 1. The formula is a string of the form "y ~ x1 / x2 / x3 / ... / xn"
/// describing the model to be fitted.
///
/// 2. y variable: Required, only 1 variable
/// Defines the distance matrix to be used, i.e. the Individuals
///
///
/// 3. x variable(s): Required, minimum 1 variable
/// Defines the grouping variables in descending order of
/// hierarchy, i.e. x1 is the highest level, x2 is the second highest level, etc.
/// The x variables can be specified as a single variable (e.g. "y ~ x") or
/// as an ordered list of variables (e.g. "y ~ x1 / x2 / x3"). Therefore the
/// minimum number of x variables required is 1.
///
/// 4. The formula must therefore contain 1 tilde and 0 or more / characters.
///

void print_help_formula(FILE* fp) {
    fprintf(fp, "Formula specification:\n");
    fprintf(fp, "  -f <formula>		: formula for AMOVA\n");
    fprintf(fp, "Formulas must have the following format:\n");
    fprintf(fp, "  <token> ~ <token> / <token> ... \n");
    fprintf(fp, "Example:\n");
    fprintf(fp, "  -f 'Individual ~ Region / Population'\n");
}

