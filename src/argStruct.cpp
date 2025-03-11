#include "shared.h"
#include "argStruct.h"
#include "io.h"

#include <time.h>
#include <htslib/hts.h> // hts_version()


//TODO ADD VCFGL style arg reading


uint8_t PROGRAM_VERBOSITY_LEVEL = 0;
char* PROGRAM_VERSION_INFO = NULL;
char* PROGRAM_COMMAND = NULL;
argStruct* args = NULL;


/// @brief usage - print usage
/// @param fp pointer to the file to print to
void help_page(FILE* fp) {

    fprintf(fp, "\n");
    fprintf(fp, "Program: ngsAMOVA\n");
    fprintf(fp, "License: GNU GPLv3.0\n");
    fprintf(fp, "Version: %s (htslib: %s)\n", NGSAMOVA_VERSION, hts_version());
    fprintf(fp, "Build: %s %s\n", __DATE__, __TIME__);
    fprintf(fp, "\n");

    fprintf(fp, "Usage: ngsAMOVA --in-(input_file_type) <input> [options]\n");

    fprintf(fp, "\n");
    fprintf(fp, "    -h, --help _____________________  Print this help message and exit\n");
    fprintf(fp, "    -v, --version __________________  Print version and build information and exit\n");
    fprintf(fp, "\n");

    fprintf(fp, "Option descriptions:\n");
    fprintf(fp, "     -s, --long-option TYPE [X] _____ Description\n");
    fprintf(fp, "     -s                               Short option (if any)\n");
    fprintf(fp, "         --long-option                Long option\n");

    fprintf(fp, "                       TYPE           Type of the argument value, can be:\n");
    fprintf(fp, "                                        - INT (integer)\n");
    fprintf(fp, "                                        - INT+ (additive integer: sum values to use together\n");
    fprintf(fp, "                                        - FLOAT (floating point number)\n");
    fprintf(fp, "                                        - STRING (string)\n");
    fprintf(fp, "                                        - FILE (filename)\n");
    fprintf(fp, "                                        - x|y|z (one of the listed values x, y or z)\n");
    fprintf(fp, "                            [X]       Default argument value (if any)\n");
    fprintf(fp, "                                _____ Connector to the option description for better readability\n");
    fprintf(fp, "\n");

    fprintf(fp, "\n");
    fprintf(fp, "General options:\n");
    fprintf(fp, "    -V, --verbose INT [0] ___________ Verbosity level\n");
    fprintf(fp, "    -@, --threads INT [1] ___________ Number of threads\n");
    fprintf(fp, "    -s, --seed INT [time] ___________ Random seed for initializing the random number generator\n");

    fprintf(fp, "\n");
    fprintf(fp, "Input/Output:\n");
    fprintf(fp, "    -i, --in-vcf FILE _______________ Input VCF/BCF file\n");
    fprintf(fp, "        --in-dm FILE ________________ Input distance matrix file\n");
    fprintf(fp, "        --in-mtd FILE _______________ Input metadata file\n");
    fprintf(fp, "        --in-jgtm FILE ______________ Input joint genotype matrix file\n");
    //TODO add in_majorminor in_blocks_tab in_blocks_bed in_regions_tab in_regions_bed


    fprintf(fp, "    -o, --output STRING ['output'] __ Output filename prefix\n");


    fprintf(fp, "\n");
    fprintf(fp, "        --print-jgtm [0]|1 __________ Print joint genotype matrix\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print joint genotype matrix (<output_prefix>.jgtm.txt)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-dm INT+ [0] _________ Print distance matrix\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print condensed distance matrix (<output_prefix>.distance_matrix.txt)\n");
    fprintf(fp, "                                   |_ 2: Print verbose distance matrix (<output_prefix>.distance_matrix.verbose.txt)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-dm-pruned INT+ [0] __ Print pruned distance matrix (requires --prune-dm 1)\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print pruned distance matrix (<output_prefix>.distance_matrix.pruned.txt)\n");
    fprintf(fp, "                                   |_ 2: Print verbose pruned distance matrix (<output_prefix>.distance_matrix.pruned.verbose.txt)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-amova INT+ [0] ______ Print AMOVA results\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print AMOVA results as CSV (<output_prefix>.amova.csv)\n");
    fprintf(fp, "                                   |_ 2: Print AMOVA results as table (<output_prefix>.amova_table.txt)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-blocks [0]|1 ________ Print tab-delimited blocks file defining the start and end positions of each block\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print blocks (<output_prefix>.blocks.txt)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-bootstrap [0]|1 _____ Print bootstrap samples\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print bootstrap samples (<output_prefix>.bootstrap_samples.tsv)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-tree [0]|1 __________ Print phylogenetic tree in Newick format\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print phylogenetic tree (<output_prefix>.newick)\n");

    fprintf(fp, "\n");
    fprintf(fp, "        --print-dxy [0]|1 ___________ Print Dxy results\n");
    fprintf(fp, "                                   |_ 0: Do not print\n");
    fprintf(fp, "                                   |_ 1: Print Dxy results (<output_prefix>.dxy.csv)\n");


    fprintf(fp, "\n");
    fprintf(fp, "Specify output file compression types:\n");
    fprintf(fp, "                                  |__ 0: No compression\n");
    fprintf(fp, "                                  |__ 1: Gzip compression (<output_prefix>.<output_suffix>.gz)\n");
    fprintf(fp, "                                  |__ 2: Bgzip compression (<output_prefix>.<output_suffix>.bgz)\n");
    fprintf(fp, "        --print-jgtm-ctype [0] ______ Print joint genotype matrix in specified compression type\n");
    fprintf(fp, "        --print-dm-ctype [0] ________ Print distance matrix in specified compression type\n");
    fprintf(fp, "        --print-dm-pruned-ctype [0] _ Print pruned distance matrix in specified compression type\n");
    fprintf(fp, "        --print-amova-ctype [0] _____ Print AMOVA results in specified compression type\n");
    fprintf(fp, "        --print-blocks-ctype [0] ____ Print blocks in specified compression type\n");


    fprintf(fp, "\n");
    fprintf(fp, "    -f, --formula STRING ____________ Formula for AMOVA\n");

    //} else if (strcasecmp("--prune-dm", arv) == 0) {
    fprintf(fp, "\n");
    fprintf(fp, "        --prune-dm [0]|1 ____________ Prune distance matrix\n");
    fprintf(fp, "                                   |_ 0: Do not prune\n");
    fprintf(fp, "                                   |_ 1: Prune distance matrix by removing individuals with the most missing data\n");

    //} else if (strcasecmp("--amova-euclid", arv)==0){
    //    // if set to 1, check distance matrix for euclidean property
    //    // if not euclidean, transform it before running AMOVA
    //    args->amova_euclid = atoi(val);
    //}
    fprintf(fp, "\n");
    fprintf(fp, "        --amova-euclid 0|[1] ________ Require Euclidean distance matrix for AMOVA\n");
    fprintf(fp, "                                   |_ 0: Do not require Euclidean distance matrix\n");
    fprintf(fp, "                                   |_ 1: Check distance matrix for Euclidean property and transform if not Euclidean for AMOVA analysis\n");

    fprintf(fp, "\n");

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

void help_page_formula(FILE* fp) {
    fprintf(fp, "Formula specification:\n");
    fprintf(fp, "  -f <formula>		: formula for AMOVA\n");
    fprintf(fp, "Formulas must have the following format:\n");
    fprintf(fp, "  <token> ~ <token> / <token> ... \n");
    fprintf(fp, "Example:\n");
    fprintf(fp, "  -f 'Individual ~ Region / Population'\n");
}


static void version_page() {
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

    return;
}


// TODO:
// check multiple of same argument
// check if given files exist here
argStruct* argStruct_get(int argc, char** argv) {


    if (argc == 0) {
        help_page(stdout);
        exit(0);
    }

    args = new argStruct;
    ASSERT(asprintf(&PROGRAM_COMMAND, "ngsAMOVA") != -1);

    while (*argv) {

        char* arv = *argv;
        char* val = *(++argv);

        if (val == NULL) {
            if (strcasecmp("--version", arv) == 0) {
                version_page();
                exit(0);
            }
            help_page(stdout);
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
//                 requires: dmat, metadata, args->formula

// -doEM <INT>    : perform EM optimization
//                 requires: gldata
//                 returns: jgtmat

// -doDxy <INT>   : estimate Dxy
//                 requires: dmat

// -doPhylo <INT> : args->do neighbor-joining tree
//                 requires: dmat

// -doMajorMinor <INT> : get major and minor alleles for each site 
//                 requires: vcfdata, alleles input file (optional)
//                 returns: alleles



        if (strcasecmp("-doUnitTests", arv) == 0) {
            args->doUnitTests = atoi(val);
            return(args);
        } else if (strcasecmp("--version", arv) == 0) {
            version_page();
            exit(0);
        } else if (strcasecmp("-doAMOVA", arv) == 0) {
            args->doAMOVA = atoi(val);
        } else if (strcasecmp("-doJGTM", arv) == 0) {
            args->doJGTM = atoi(val);
        } else if (strcasecmp("-doEM", arv) == 0) {
            args->doEM = atoi(val);
            // when performing args->block bootstrapping em, use --verbose >1 to print per em optimization termination reason, last d value, and number of iterations
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
        } else if (strcasecmp("-doDryRun", arv) == 0 || strcasecmp("--dryrun", arv) == 0) {
            args->doDryRun = atoi(val);
            //TODO implement
        } else if (strcasecmp("--alloc", arv) == 0) {
            args->alloc_strategy = atoi(val);
            //TODO RMME?
        } else if (strcasecmp("--bcf-src", arv) == 0) {
            args->bcfSrc = atoi(val);
        } else if (strcasecmp("--prune-dm", arv) == 0) {
            // any args->downstream analysis will use the pruned version of the dmat
            // does not effect any printing, its printing is controlled by:
            // --print-pruned-dm 1 will print pruned dm
            // see ARG_PRINT_PRUNED_DM
            args->prune_dmat = atoi(val);
        } else if (strcasecmp("--amova-euclid", arv) == 0) {
            // if set to 1, check distance matrix for euclidean property
            // if not euclidean, transform it before running AMOVA
            args->amova_euclid = atoi(val);
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

        else if ((strcasecmp("--in-vcf", arv) == 0) || (strcasecmp("-i", arv) == 0)) {
            args->in_vcf_fn = strdup(val);
        } else if ((strcasecmp("--in-dm", arv) == 0)) {
            args->in_dm_fn = strdup(val);
        } else if ((strcasecmp("--in-jgtm", arv) == 0)) {
            args->in_jgtmat_fn = strdup(val);
        } else if ((strcasecmp("--in-majorminor", arv) == 0)) {
            args->in_majorminor_fn = strdup(val);
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
        // #################################################################
        //
        // Use printing commands to specify the output files to be generated.
        // This is only needed for output files that are not the default output files
        //   associated with the analyses specified by the action commands.
        //
        // -> --print-jgtm
        else if (strcasecmp("--print-jgtm", arv) == 0) {
            args->print_jgtm = atoi(val);
        } else if (strcasecmp("--print-jgtm-ctype", arv) == 0) {
            args->print_jgtm_ctype = atoi(val);

            // -> --print-dm
        } else if (strcasecmp("--print-dm", arv) == 0) {
            args->print_dm = atoi(val);
        } else if (strcasecmp("--print-dm-ctype", arv) == 0) {
            args->print_dm_ctype = atoi(val);

            // -> --print-pruned-dm
        } else if (strcasecmp("--print-pruned-dm", arv) == 0) {
            args->print_pruned_dm = atoi(val);
        } else if (strcasecmp("--print-pruned-dm-ctype", arv) == 0) {
            args->print_pruned_dm_ctype = atoi(val);

            // -> --print-amova
        } else if (strcasecmp("--print-amova", arv) == 0) {
            args->print_amova = atoi(val);
        } else if (strcasecmp("--print-amova-ctype", arv) == 0) {
            args->print_amova_ctype = atoi(val);

            // -> --print-blocks
        } else if (strcasecmp("--print-blocks", arv) == 0) {
            args->print_blocks = atoi(val);
        } else if (strcasecmp("--print-blocks-ctype", arv) == 0) {
            args->print_blocks_ctype = atoi(val);



            // -> --print-bootstrap
                // prints to <prefix>.bootstrap_samples.tsv
                // with header:
                // Rep\tPos\targs->blockID\targs->blockContig\targs->blockStart\targs->blockEnd
                // Rep 	Pos 	args->blockID 	args->blockContig args->blockStart 	args->blockEnd
                // Rep:        Replicate number
                // Pos:        Position of the sampled args->block in the replicates synthetic genome (0-based)
                // args->blockID:    ID of the sampled args->block
                // args->blockContig: Name of the chromosome to which the args->block belongs to
                // args->blockStart: 1-based, inclusive [start position of the args->block with the given args->blockID
                // args->blockEnd:   1-based, inclusive end] position of the args->block with the given args->blockID
                //
                // e.g.
                // if we have 4 args->blocks in the original genome, each with a size of 1000 bps 
                // {args->block1, args->block2, args->block3, args->block4} 
                // (0     , 1     , 2     , 3) args->blockIDs
                // internal representation start end pos: (0-based inclusive start exclusive end)
                // (all belonging to chr1)
                // args->block1 start:0,    end:1000 
                // args->block2 start:1000, end:2000
                // args->block3 start:2000, end:3000
                // args->block4 start:3000, end:4000
                // args->blocks tab representation: (1-based inclusive start end)
                // args->block1 start:1,    end:1000
                // args->block2 start:1001, end:2000
                // args->block3 start:2001, end:3000
                // args->block4 start:3001, end:4000
                //
                // and we have one replicate which has sampled the following ordered set:
                // {args->block3, args->block1, args->block2, args->block4}
                //
                // Our output file output.bootstrap_samples.tsv will have:
                // Rep 	Pos 	args->blockID 	args->blockContig args->blockStart 	args->blockEnd
                // 0    0       2           chr1        2001        3000
                // 0    1       0           chr1        0           1000
                // 0    2       1           chr1        1001        2000
                // 0    3       3           chr1        3001        4000
                // 

        } else if (strcasecmp("--print-bootstrap", arv) == 0) {
            args->print_bootstrap = atoi(val);
        } else if (strcasecmp("--print-bootstrap-ctype", arv) == 0) {
            args->print_bootstrap_ctype = atoi(val);
        } else if (strcasecmp("--print-tree", arv) == 0) {
            args->print_tree = atoi(val);
        } else if (strcasecmp("--print-tree-ctype", arv) == 0) {
            args->print_tree_ctype = atoi(val);
        } else if (strcasecmp("--print-dxy", arv) == 0) {
            args->print_dxy = atoi(val);
        } else if (strcasecmp("--print-dxy-ctype", arv) == 0) {
            args->print_dxy_ctype = atoi(val);
        } else if (strcasecmp("--print-ibd", arv) == 0) {
            args->print_ibd = atoi(val);
        } else if (strcasecmp("--print-ibd-ctype", arv) == 0) {
            args->print_ibd_ctype = atoi(val);
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
        // #    args->block BOOTSTRAPPING COMMANDS                                 #
        // #    [--block-size/-bs] [--blocks-tab] [--blocks-bed]             #
        // ###################################################################
        //
        // Use args->block bootstrapping commands to specify the args->blocks to be used in the args->block bootstrapping analyses.
        //
        // There are two types of args->block bootstrapping specifications:
        // - args->block size specification
        // - args->block list specification
        //
        // args->block size specification commands are of the form `--block-size <int>` or `-bs <int>`
        //  where <int> is the size of the args->blocks. Using the VCF file as input, the args->blocks are enumerated
        //  by reading the contig sizes from the VCF header and dividing the contigs into args->blocks of size <int>.
        //
        // args->block list specification commands are of the form `--blocks-tab <filename>` or `-bf <filename>`
        //  where <filename> is the name of the file containing the list of args->blocks.
        //
        // The args->blocks file should be in one of the following formats:
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
        // -> Only one args->block bootstrapping specification command can be used at a time.

        else if (strcasecmp("--block-size", arv) == 0 || strcasecmp("-bs", arv) == 0) {
            // read args->block size as float and convert to int
            // to allow for the use of scientific notation (e.g. 1e6)
            args->blockSize = (int)atof(val);
        }

        else if (strcasecmp("--blocks-tab", arv) == 0) {
            args->in_blocks_tab_fn = strdup(val);
        } else if (strcasecmp("--blocks-bed", arv) == 0) {
            args->in_blocks_bed_fn = strdup(val);
        } else if (strcasecmp("--a2f", arv) == 0) {
            args->in_a2f_fn = strdup(val);


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
        else if (strcasecmp("--min-a2f", arv) == 0) {
            args->min_info_a2f = atof(val);
        }

        else if (strcasecmp("--min-a2c", arv) == 0) {
            args->min_a2c = atoi(val);
            if (args->min_a2c < 2) {
                ERROR("--min-a2c is set to %d. Allowed range of values is 2<value<N", args->min_a2c);
            }
        }

        else if (strcasecmp("--errorprop", arv) == 0) {
            args->ibd_errorprop = atof(val);
        }

        else if (strcasecmp("--errormax", arv) == 0) {
            args->ibd_errormax = atof(val);
        }

        else if (strcasecmp("--ibd-alpha", arv) == 0) {
            args->ibd_alpha = atof(val);
            //TODO add value check

        } else if (strcasecmp("--ibd-beta", arv) == 0) {
            args->ibd_beta = atof(val);
            //TODO add value check

        } else if (strcasecmp("--ibd-dynamic-alpha", arv) == 0) {
            args->ibd_dynamic_alpha = atoi(val);
            //TODO add value check

        } else if (strcasecmp("--ibd-dynamic-beta", arv) == 0) {
            args->ibd_dynamic_beta = atoi(val);
            //TODO add value check

        } else if (strcasecmp("--ibdlod", arv) == 0) {
            args->ibd_ibdlod = atof(val);
            if (args->ibd_ibdlod < 0) {
                ERROR("--ibdlod must be a positive number.");
            }
        } else if (strcasecmp("--ibdtrim", arv) == 0) {
            args->ibd_ibdtrim = atof(val);
            if (args->ibd_ibdtrim < 0) {
                ERROR("--ibdtrim must be a non-negative number.");
            }

            // add opt: max-missing
        } else if (strcasecmp("--max-missing", arv) == 0) {
            NEVER;
            //TODO implement
            int tmp = atoi(val);
            if (tmp < 0) {
                ERROR("--max-missing must be a non-negative number.");
            }
            args->ibd_segment_max_n_missing_sites = tmp;
        } else if (strcasecmp("--em-tole", arv) == 0)
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

        } else if (strcasecmp("--rm-multiallelic-sites", arv) == 0) {
            args->rmMultiallelicSites = atoi(val);

        } else if (strcasecmp("--allow-mispairs", arv) == 0) {
            // if an individual pair has no shared sites, should the program drop the pair or give err and exit?
            args->allow_mispairs = atoi(val);

        } else if (strcasecmp("--min-pairsites", arv) == 0) {
            // drop the pair if the number of shared sites is less than this value
            args->pair_min_n_sites = (int)atof(val);

        } else if (strcasecmp("--min-npairs", arv) == 0) {
            args->min_n_pairs = atoi(val);

        } else if (strcasecmp("-h", arv) == 0 || strcasecmp("--help", arv) == 0) {
            help_page(stdout);
            exit(0);

        } else {
            ERROR("Unknown argument: \'%s\'\n", arv);
        }

        char* temp = NULL;
        ASSERT(asprintf(&temp, "%s %s %s", PROGRAM_COMMAND, arv, val) != -1);
        if (PROGRAM_COMMAND != NULL) {
            FREE(PROGRAM_COMMAND);
        }
        PROGRAM_COMMAND = temp;

        ++argv;
    }

    LOG("Running command: %s", PROGRAM_COMMAND);

    if (args->doUnitTests) {
        return(args);
    }

    //--------------------------------------------
    // -> check printing options
    // check for unset printing options and set defaults

    // print jgtm
    if (args->print_jgtm == -1) {
        // default: do not print jgtm (ARG_PRINT_NOPRINT)
        args->print_jgtm = ARG_PRINT_NOPRINT;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_jgtm, ARG_PRINT_NOPRINT, ARG_PRINT_JGTM, "--print-jgtm");
    }
    if (args->print_jgtm_ctype == -1) {
        // default: print jgtm without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_jgtm_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_jgtm_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-jgtm-ctype");
    }

    // print dm
    if (args->print_dm == -1) {
        // default: do not print dm (ARG_INTPLUS_UNSET)
        args->print_dm = ARG_INTPLUS_UNSET;
    } else {
        CHECK_ARG_INTERVAL_INTPLUS(args->print_dm, ARG_INTPLUS_PRINT_DM_VERBOSE, "--print-dm");
    }
    if (args->print_dm_ctype == -1) {
        // default: print dm without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_dm_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_dm_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-dm-ctype");
    }

    // print pruned dm
    if (args->print_pruned_dm == -1) {
        // default: do not print pruned dm (ARG_INTPLUS_UNSET)
        args->print_pruned_dm = ARG_INTPLUS_UNSET;
    } else {
        CHECK_ARG_INTERVAL_INTPLUS(args->print_pruned_dm, ARG_INTPLUS_PRINT_PRUNED_DM_VERBOSE, "--print-pruned-dm");
    }
    if (args->print_pruned_dm_ctype == -1) {
        // default: print pruned dm without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_pruned_dm_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_pruned_dm_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-pruned-dm-ctype");
    }

    // print amova
    if (args->print_amova == -1) {
        // default: do not print amova (ARG_INTPLUS_UNSET)
        args->print_amova = ARG_INTPLUS_UNSET;
    } else {
        CHECK_ARG_INTERVAL_INTPLUS(args->print_amova, ARG_INTPLUS_PRINT_AMOVA_TABLE, "--print-amova");
    }
    if (args->print_amova_ctype == -1) {
        // default: print amova without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_amova_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_amova_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-amova-ctype");
    }

    // print blocks
    if (args->print_blocks == -1) {
        // default: do not print args->blocks (ARG_PRINT_NOPRINT)
        args->print_blocks = ARG_PRINT_NOPRINT;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_blocks, ARG_PRINT_NOPRINT, ARG_PRINT_BLOCKS, "--print-blocks");
    }
    if (args->print_blocks_ctype == -1) {
        // default: print args->blocks without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_blocks_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_blocks_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-blocks-ctype");
    }

    // print bootstrap samples
    if (args->print_bootstrap == -1) {
        // default: do not print bootstrap samples (ARG_PRINT_NOPRINT)
        args->print_bootstrap = ARG_PRINT_NOPRINT;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_bootstrap, ARG_PRINT_NOPRINT, ARG_PRINT_BOOTSTRAP, "--print-bootstrap");
    }
    if (args->print_bootstrap_ctype == -1) {
        // default: print bootstrap samples without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_bootstrap_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_bootstrap_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-bootstrap-ctype");
    }

    // print tree
    if (args->print_tree == -1) {
        // default: do not print tree (ARG_PRINT_NOPRINT)
        args->print_tree = ARG_PRINT_NOPRINT;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_tree, ARG_PRINT_NOPRINT, ARG_PRINT_TREE, "--print-tree");
    }
    if (args->print_tree_ctype == -1) {
        // default: print tree without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_tree_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_tree_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-tree-ctype");
    }

    // print dxy
    if (args->print_dxy == -1) {
        // default: do not print dxy (ARG_PRINT_NOPRINT)
        args->print_dxy = ARG_PRINT_NOPRINT;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_dxy, ARG_PRINT_NOPRINT, ARG_PRINT_DXY, "--print-dxy");
    }
    if (args->print_dxy_ctype == -1) {
        // default: print dxy without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_dxy_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_dxy_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-dxy-ctype");
    }

    // print ibd segments
    if (args->print_ibd == -1) {
        // default: do not print ibd (ARG_INTPLUS_UNSET)
        args->print_ibd = ARG_INTPLUS_UNSET;
    } else {
        CHECK_ARG_INTERVAL_INTPLUS(args->print_ibd, ARG_INTPLUS_PRINT_IBD_PERSITE_SMOOTHED_IBD_SCORES, "--print-ibd");
    }
    if (args->print_ibd_ctype == -1) {
        // default: print ibd without compression (PROGRAM_OUTFILE_CTYPE_RAW)
        args->print_ibd_ctype = PROGRAM_OUTFILE_CTYPE_RAW;
    } else {
        CHECK_ARG_INTERVAL_II_INT(args->print_ibd_ctype, PROGRAM_OUTFILE_CTYPE_RAW, PROGRAM_OUTFILE_CTYPE_BGZ, "--print-ibd-ctype");
    }
    //--------------------------------------------

    if (NULL == args->out_fnp) {
        args->out_fnp = strdup("output");
        fprintf(stderr, "\n\t-> -out <output_prefix> not set; will use %s as a prefix for output files.\n", args->out_fnp);
    } else {
        fprintf(stderr, "\n\t-> -out <output_prefix> is set to %s. Output files will have this prefix.\n", args->out_fnp);
    }


    if (-1 != args->allow_mispairs) {
        CHECK_ARG_INTERVAL_01(args->allow_mispairs, "--allow-mispairs");
        WARN("(--allow-mispairs %d) For downstream analyses, use this option with caution as it may lead to decreased power.", args->allow_mispairs);
    } else {
        args->allow_mispairs = 0;
        LOG("--allow-mispairs is not set, setting to default value %d (do not drop pairs). Program will give error and exit if an individual pair has no shared sites.", args->allow_mispairs);
    }

    if (-1 != args->pair_min_n_sites) {
        CHECK_ARG_INTERVAL_INT(args->pair_min_n_sites, 1, INT_MAX, "--min-pairsites");
    } else {
        args->pair_min_n_sites = 1;
        LOG("--min-pairsites is not set, setting to default value %d (perform the action defined via --allow-mispairs if an individual pair has 0 shared sites).", args->pair_min_n_sites);
    }

    if (args->print_pruned_dm) {
        if (args->prune_dmat != 1) {
            ERROR("Cannot print pruned distance matrix when pruning is not enabled.");
        }
    }


    if (-1 != args->min_n_pairs) {
        // requires: --allow-mispairs 1
        CHECK_ARG_INTERVAL_INT(args->min_n_pairs, 0, INT_MAX, "--min-npairs");
        // if (0 == args->allow_mispairs) {
        //     ERROR("Cannot set --min-npairs without setting --allow-mispairs to 1. Please set --allow-mispairs to 1.");
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
        LOG("Input is VCF file: %s", args->in_vcf_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_VCF;
        if (0 == args->bcfSrc) {
            ERROR("BCF data source is necessary for VCF input. Please set the BCF data source using --bcf-src.");
        }
    }

    if (NULL != args->in_dm_fn) {
        LOG("Input is distance matrix file: %s", args->in_dm_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_DM;
    }

    if (NULL != args->in_mtd_fn) {
        //fprintf(stderr, "\n[INFO]\tFound input metadata file: %s\n", args->in_mtd_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_METADATA;
    }

    if (NULL != args->in_majorminor_fn) {
        //fprintf(stderr, "\n[INFO]\tFound input major/minor alleles file: %s\n", args->in_majorminor_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_MAJORMINOR;
    }

    if (NULL != args->in_blocks_bed_fn) {
        //fprintf(stderr, "\n[INFO]\tFound input args->blocks BED file: %s\n", args->in_blocks_bed_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_BLOCKS;
    }

    if (NULL != args->in_blocks_tab_fn) {
        //fprintf(stderr, "\n[INFO]\tFound input args->blocks TSV file: %s\n", args->in_blocks_tab_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_BLOCKS;
    }

    if (NULL != args->in_regions_bed_fn) {
        //fprintf(stderr, "\n[INFO]\tFound input regions BED file: %s\n", args->in_regions_bed_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_REGIONS;
    }

    if (NULL != args->in_regions_tab_fn) {
        //TODO test
        //fprintf(stderr, "\n[INFO]\tFound input regions TSV file: %s\n", args->in_regions_tab_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_REGIONS;
    }

    if (NULL != args->in_jgtmat_fn) {
        //fprintf(stderr, "\n[INFO]\tFound input joint genotype matrix file: %s\n", args->in_jgtmat_fn);
        LOG("Input is joint genotype matrix file: %s", args->in_jgtmat_fn);
        args->in_ft = args->in_ft | ARG_INTPLUS_INPUT_JGTMAT;
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

    if (args->print_jgtm > 0) {
        if (!PROGRAM_HAS_INPUT_VCF) {
            ERROR("Cannot print joint genotype matrix without VCF input.");
        }
    }

    if (PROGRAM_WILL_USE_BCF_FMT_GT) {
        if (args->doDist) {
            if (!args->doJGTM) {
                ERROR("Cannot calculate distance matrix without calculating joint genotype matrix.");
            }
        }
    }

    if (PROGRAM_HAS_INPUT_METADATA && (!PROGRAM_NEEDS_METADATA)) {
        WARN("Metadata file is provided but no analysis requires it; will ignore the metadata file %s.", args->in_mtd_fn);
        if (NULL != args->formula) {
            WARN("args->formula is set but no analysis requires it; will ignore the args->formula '%s'.", args->formula);
        }
    }

    if ((!(PROGRAM_HAS_INPUT_METADATA)) && PROGRAM_NEEDS_METADATA) {
        ERROR("Metadata file is required for the specified analyses. Please provide the metadata file using --metadata.");
    }

    if (PROGRAM_HAS_INPUT_METADATA && PROGRAM_NEEDS_METADATA) {
        if (args->formula == NULL) {
            ERROR("Metadata file is provided but no args->formula is set. Please set the args->formula using -f/--formula.");
        }
    }


    if (args->doBlockBootstrap) {
        if (PROGRAM_WILL_USE_BCF_FMT_GT) {
            ERROR("Program cannot use GT data for block bootstrapping.");
        }


        if (args->nBootstraps != 0) {
            CHECK_ARG_INTERVAL_INT(args->nBootstraps, 1, 100000, "-nb/--nBootstraps");
        }else{
            ERROR("-doBlockBootstrap is set but the number of bootstraps is not set. Please set the number of bootstraps using --nBootstraps.");
        }

        if (-1.0 == args->bootstrap_pctci) {
            args->bootstrap_pctci = 95.0;
            LOG("Bootstrap confidence interval is not set, setting to default value %.17g %%.", args->bootstrap_pctci);

        } else {
            CHECK_ARG_INTERVAL_EE_DBL(args->bootstrap_pctci, 0.0, 100.0, "--bootstrap-ci");
        }

        if (args->blockSize == -1 && args->in_blocks_tab_fn == NULL && args->in_blocks_bed_fn == NULL) {
            ERROR("Block bootstrapping is enabled but no block specification is provided. Please provide block size (--block-size) or a blocks specification file (--in-blocks-tab or --in-blocks-bed).");
        }
        // TODO in the new arg checking system, directly check if unset instead on relying on the unset values?
        /// like if getArg(args, "--nBootstraps", &args->nBootstraps) returns false, then it is unset
        // if (args->nBootstraps == 0) {
        //     ERROR("args->block bootstrapping is enabled but the number of bootstraps is not set. Please set the number of bootstraps using -nb/--nBootstraps.");
        // }

    }

    if (PROGRAM_HAS_INPUT_METADATA && (!PROGRAM_NEEDS_METADATA)) {
        ERROR("Metadata file is provided but none of the specified analyses use it. Please remove the metadata file or set the correct analysis.");

    }

    ASSERT(asprintf(&PROGRAM_VERSION_INFO, "ngsAMOVA [version: %s] [build: %s %s] [htslib: %s]", NGSAMOVA_VERSION, __DATE__, __TIME__, hts_version()) != -1);


    if (args->doIbd == ARG_DOIBD_GT_METHOD) {
        if (args->min_a2c < 2) {
            ERROR("%s requires %s to be at least 2.", "--doIbd", "--min-a2c");
        }
    }

    if (args->doIbd) {

        const double e = args->ibd_errormax;
        const double x = 1.0 - e;

        args->ibd_max_error_array[0] = x * x * x * x;
        args->ibd_max_error_array[1] = (e) * (x * x * x);
        args->ibd_max_error_array[2] = (e * e) * (x * x);
        args->ibd_max_error_array[3] = (e * e * e) * (x);
        args->ibd_max_error_array[4] = (e * e * e * e);
    }

    if (args->minInd == -1) {
        // LOG("args->minInd is not set. Default is setting args->minInd to 2; will use sites that is nonmissing for both individuals in each individual pair.");
        // minInd filter is not set
        // args->minInd = 2;
    } else if (args->minInd == 2) {
        LOG("--minInd is set to 2; will use sites that is nonmissing for both individuals in a pair.");
    } else if (args->minInd == 1 || args->minInd < -1) {
        ERROR("--minInd is set to %d. Minimum value allowed for --minInd is 2.", args->minInd);
    }

    if (args->seed == -1) {
        args->seed = time(NULL);
        srand48(args->seed);
        LOG("--seed is not defined, will use current time as random seed for the random number generator. seed is now set to: %d.\n", args->seed);
        WARN("Used the current time as random seed for random number generator. For parallel runs this may cause seed collisions. Hence, it is recommended to set the seed manually using `--seed <INTEGER>`.");
    } else {
        srand48(args->seed);
        LOG("--seed is set to: %d.", args->seed);
    }

    if (args->in_dm_fn != NULL) {
        if (args->doEM != 0) {
            fprintf(stderr, "\n[ERROR]\t-doEM %i cannot be used with -in_dm %s.", args->doEM, args->in_dm_fn);
            exit(1);
        }
    }

    if (args->windowSize != 0) {
        fprintf(stderr, "\n[INFO]\t-> -windowSize %d; will use sliding winargs->dows of size %d\n", args->windowSize, args->windowSize);
        NEVER;
    }

    //----------------------------------------------------------------------------------//
    // -doDist      defines the method to estimate the pairwise distance matrix
    // (get dmat)
    //              default: 0
    //
    //              0: args->do not estimate distance matrix
    //              1: calculate the pairwise distance matrix using the method defined via --dm-method
    //              2: read the distance matrix from the file defined via --in-dm
    //              3: args->do args->doDist 1 for original dataset and the args->block bootstrapped datasets

    // check for unset values and set defaults
    if (args->dm_transform == -1) {
        args->dm_transform = DMAT_INTPLUS_TRANSFORM_NONE;
        LOG("Distance matrix transform is not set, setting to default value %d (none).", args->dm_transform);
    } else {
        CHECK_ARG_INTERVAL_INTPLUS(args->dm_transform, DMAT_INTPLUS_TRANSFORM_MAX, "--dm-transform");
    }

    if (-1 == args->dm_method) {
        args->dm_method = DMAT_METHOD_DIJ;
        LOG("Distance matrix method is not set, setting to default value %d (Dij).", args->dm_method);
    } else {
        CHECK_ARG_INTERVAL_INT(args->dm_method, 0, 9, "--dm-method");
    }

    if (0 == args->doDist) {
        ;;
    } else if ((1 == args->doDist) || (3 == args->doDist)) {
        ;;

    } else if (2 == args->doDist) {
        if (NULL == args->in_dm_fn) {
            ERROR("-doDist %d requires a distance matrix file (--in-dm <file>).", args->doDist);
        }
        LOG("-doDist %d; will read the distance matrix from the file %s.", args->doDist, args->in_dm_fn);
    } else {
        ERROR("-doDist %d is not a valid option.", args->doDist);
    }


    //----------------------------------------------------------------------------------//
    // -doAMOVA
    //             default: 0
    //             0: args->do not perform AMOVA
    //             1: perform AMOVA 
    //             2: perform AMOVA with args->block bootstrapping (requires: args->block size, args->block file etc)
    //             3: perform AMOVA with permutation test (orig framework, requires: npermut)
    // requires:
    // - args->formula
    // - metadata
    // - distance matrix

    if (args->doAMOVA) {
        if (ARG_DOAMOVA_SINGLERUN == args->doAMOVA) {
            IO::requireArgStr(args->formula, "--formula/-f", "-doAMOVA 1");
            IO::requireArgFile(args->in_mtd_fn, "--metadata/-m", "-doAMOVA 1");

            if (0 == args->doDist && NULL == args->in_dm_fn) {
                ERROR("-doAMOVA %d requires either (1) a method to calculate the pairwise distance matrix (-doDist <int>) or (2) a distance matrix file (--in-dm <file>).", args->doAMOVA);
            }

        } else if (ARG_DOAMOVA_BOOTSTRAP == args->doAMOVA) {
            IO::requireArgStr(args->formula, "--formula/-f", "-doAMOVA 1");
            IO::requireArgFile(args->in_mtd_fn, "--metadata/-m", "-doAMOVA 1");

            if (0 == args->doDist && NULL == args->in_dm_fn) {
                ERROR("-doAMOVA %d requires either (1) a method to calculate the pairwise distance matrix (-doDist <int>) or (2) a distance matrix file (--in-dm <file>).", args->doAMOVA);
            }
            if (!PROGRAM_WILL_USE_BCF_FMT_GL) {
                ERROR("-doAMOVA %d requires genotype likelihoods in the input VCF file (--bcf-src %d).", args->doAMOVA, ARG_INTPLUS_BCFSRC_FMT_GL);
            }
        } else if (ARG_DOAMOVA_PERMTEST) {
            NEVER;//TODO 
        } else {
            ERROR("-doAMOVA %d is not a valid option.", args->doAMOVA);
        }
    }

    // If the input data source is genotype likelihoods(`- - bcf - src 1`), the EM optimization(`-doEM 1`) is needed to obtain the joint genotype matrix.
    if (0 == args->doEM) {
        //
        if (-1 != args->tole) {
            ERROR("--em-tole %e requires -doEM != 0.", args->tole);
        }
        if (-1 != args->maxEmIter) {
            ERROR("--maxEmIter %d requires -doEM != 0.", args->maxEmIter);
        }
    } else if (1 == args->doEM) {
        if (-1 == args->tole) {
            args->tole = 1e-6;
            LOG("-tole is not set, setting to default value %e. Will terminate the EM algorithm if the change in the log-likelihood is less than %e.", args->tole, args->tole);
        } else {
            LOG("-tole is set to %e. Will terminate the EM algorithm if the change in the log-likelihood is less than %e.", args->tole, args->tole);
        }

        if (-1 == args->maxEmIter) {
            args->maxEmIter = 500;
            LOG("-maxEmIter is not set, setting to default value %d. Will terminate the EM algorithm if the number of iterations exceed %d.", args->maxEmIter, args->maxEmIter);
        } else {
            LOG("-maxEmIter is set to %d. Will terminate the EM algorithm if the number of iterations exceed %d.", args->maxEmIter, args->maxEmIter);
        }

    } else {
        ERROR("-doEM %d is not a valid option.", args->doEM);
    }

    //----------------------------------------------------------------------------------//
// -doDxy
    //              default: 0
    //              0: args->do not estimate dxy
    //              1: estimate dxy from the distance matrix for all groups in all hierarchical levels defined in the metadata file (requires: method to obtain distance matrix, metadata file)

    // : (args->doDist 1,2 or args->in_ft == args->in_DM), (--metadata-file <METADATA_FILE>)

    if (0 == args->doDxy) {
        //
    } else if (1 == args->doDxy) {
        IO::requireArgFile(args->in_mtd_fn, "--metadata/-m", "-doDxy 1");
        if (args->in_mtd_fn == NULL) {
            ERROR("-doDxy %d requires --metadata/-m <file>.", args->doDxy);
        }
    } else {
        ERROR("-doDxy %d is not a valid option.", args->doDxy);
    }














    //----------------------------------------------------------------------------------//
    // -doPhylo
    //              default: 0
    //              0: args->do not run phylogenetic tree construction
    //              1: construct phylogenetic tree using neighbor joining with individuals as leaf nodes
    //              2: construct phylogenetic tree using neighbor joining with groups as leaf nodes (requires: `-doDxy`)

    if (0 == args->doPhylo) {
        //
    } else if (1 == args->doPhylo) {
        // if (0 == args->doDist && NULL == args->in_dm_fn) {
        //     ERROR("-doPhylo %d requires a distance matrix (either -doDist <int> or --in-dm <file>).", args->doPhylo);
        // }
    } else if (2 == args->doPhylo) {
    } else {
        ERROR("-doPhylo %d is not a valid option.", args->doPhylo);
    }

    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GL) {
        LOG("%s INT+ argument value contains %d. Program will use the GL tag from input VCF file.", "--bcf-src", ARG_INTPLUS_BCFSRC_FMT_GL);
    }
    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GT) {
        LOG("%s INT+ argument value contains %d. Program will use the GT tag from input VCF file.", "--bcf-src", ARG_INTPLUS_BCFSRC_FMT_GT);
    }

    // -------------------------------
    // args->doMajorMinor
    if (args->in_majorminor_fn != NULL) {
        if (args->doMajorMinor != ARG_DOMAJORMINOR_INFILE) {
            ERROR("%s requires %s to be set to %d.", "--in-majorminor", "--doMajorMinor", ARG_DOMAJORMINOR_INFILE);
        }
    }
    if (args->doMajorMinor == ARG_DOMAJORMINOR_BCF_REFALT1) {
        if (args->doIbd == ARG_DOIBD_GT_METHOD) {
            LOG("%s is set to %d, and %s is set to %d. Program will attempt using the REF allele as the major allele and the ALT1 allele as the minor allele. N.B. With this option, program expects REF allele to be the allele with the highest allele count and the first ALT allele to be the allele with the second highest allele count. If this is not the case, the program will give an error and exit.", "-doMajorMinor", ARG_DOMAJORMINOR_BCF_REFALT1, "-doIbd", ARG_DOIBD_GT_METHOD);
        } else {
            LOG("%s is set to %d. Program will use the REF allele as the major allele and the the first ALT allele as the minor allele.", "-doMajorMinor", ARG_DOMAJORMINOR_BCF_REFALT1);
        }
        if (args->in_majorminor_fn != NULL) {
            ERROR("Program cannot use both -doMajorMinor %d and --in-majorminor <file>.", args->doMajorMinor);
        }

    } else if (args->doMajorMinor == ARG_DOMAJORMINOR_INFILE) {
        if (args->in_majorminor_fn == NULL) {
            ERROR("-doMajorMinor %d requires --in-majorminor <file> to be set.", args->doMajorMinor);
        }
        if (args->in_majorminor_fn != NULL) {
            LOG("%s is set to %d, and %s is provided. Program will use the major and minor alleles from the file %s.", "-doMajorMinor", ARG_DOMAJORMINOR_INFILE, "--in-majorminor", args->in_majorminor_fn);
        }
    } else if (args->doMajorMinor == ARG_DOMAJORMINOR_BCF_AC) {
        LOG("%s is set to %d. Program will use the allele counts from the BCF file to determine the major and minor alleles. N.B. This option requires the BCF file to contain the genotypes (FORMAT/GT tag) or allele counts (INFO/AC tag).", "-doMajorMinor", ARG_DOMAJORMINOR_BCF_AC);
    }

    if (ARG_DOMAJORMINOR_UNSET == args->doMajorMinor && NULL == args->in_majorminor_fn) {
        if (args->doIbd) {
            ERROR("Program cannot perform IBD args->block detection without major/minor allele information (-doMajorMinor).");
        }
        if (args->doEM) {
            NEVER;//TODO 10gl not implemented yet
        }
    }

    // -------------------------------

        // TODO using regions with args->block definitions
        // TODO handle empty args->blocks
    if (args->blockSize > 0 && args->in_blocks_tab_fn != NULL) {
        fprintf(stderr, "\n[ERROR]\t-> `--block-size` cannot be used with `--in-blocks-tab`.\n");
        exit(1);
    }
    if (args->blockSize > 0 && args->in_blocks_bed_fn != NULL) {
        fprintf(stderr, "\n[ERROR]\t-> `--block-size` cannot be used with `--in-blocks-bed`.\n");
        exit(1);
    }

    if (args->in_blocks_tab_fn != NULL || args->in_blocks_bed_fn != NULL) {
        if (args->in_regions_tab_fn != NULL || args->in_regions_bed_fn != NULL || args->in_region != NULL) {
            fprintf(stderr, "\n[ERROR]\targs->block definitions cannot be used with region definitions, yet.\n");
            exit(1);
        }
    }

    if (args->doEM) {
        if (0 == args->doJGTM) {
            ERROR("-doEM requires -doJGTM");
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



    if (args->doDryRun) {
        if (!PROGRAM_HAS_INPUT_VCF) {
            ERROR("(%s) Dry run is only available for VCF input. Exiting...", "--dryrun");
        }
    }

    return args;
}

void argStruct_destroy(argStruct* _args) {

    if (_args->in_vcf_fn != NULL) {
        FREE(_args->in_vcf_fn);
    }
    if (_args->in_dm_fn != NULL) {
        FREE(_args->in_dm_fn);
    }
    if (_args->in_mtd_fn != NULL) {
        FREE(_args->in_mtd_fn);
    }
    if (_args->in_jgtmat_fn != NULL) {
        FREE(_args->in_jgtmat_fn);
    }

    if (_args->out_fnp != NULL) {
        FREE(_args->out_fnp);
    }

    if (_args->in_region != NULL) {
        FREE(_args->in_region);
    }
    if (_args->in_regions_tab_fn != NULL) {
        FREE(_args->in_regions_tab_fn);
    }
    if (_args->in_regions_bed_fn != NULL) {
        FREE(_args->in_regions_bed_fn);
    }

    if (_args->in_blocks_tab_fn != NULL) {
        FREE(_args->in_blocks_tab_fn);
    }
    if (_args->in_blocks_bed_fn != NULL) {
        FREE(_args->in_blocks_bed_fn);
    }

    if (_args->in_majorminor_fn != NULL) {
        FREE(_args->in_majorminor_fn);
    }

    if (_args->in_a2f_fn != NULL) {
        FREE(_args->in_a2f_fn);
    }

    if (_args->command != NULL) {
        FREE(_args->command);
    }

    if (_args->formula != NULL) {
        FREE(_args->formula);
    }

    if (NULL != args->a2freqs) {
        FREE(args->a2freqs);
    }

    if (NULL != PROGRAM_VERSION_INFO) {
        FREE(PROGRAM_VERSION_INFO);
    }
    if (NULL != PROGRAM_COMMAND) {
        FREE(PROGRAM_COMMAND);
    }

    ks_free(&_args->outfiles_list);

    delete _args;
}



