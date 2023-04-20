#include "argStruct.h"
#include "paramStruct.h"
#include "io.h"

// default: 0 (verbose mode off)
u_char VERBOSE = 0;

// TODO check multiple of same argument
/// @brief argStruct_get read command line arguments
/// @param argc
/// @param argv
/// @return pointer to argStruct
argStruct *argStruct_get(int argc, char **argv)
{

	argStruct *args = new argStruct;

	while (*argv)
	{

		// TODO check if given files exist here

		char *arv = *argv;
		char *val = *(++argv);

		if (val == NULL)
		{
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
		//   -doNJ <int> : do neighbor-joining
		//   -doDist <int> : calculate pairwise distances

		if (strcasecmp("-doAMOVA", arv) == 0)
		{
			args->doAMOVA = atoi(val);
		}
		else if (strcasecmp("-doEM", arv) == 0)
		{
			args->doEM = atoi(val);
		}

		else if (strcasecmp("-doDxy", arv) == 0)
		{
			if (strIsNumeric(val))
			{
				args->doDxy = atoi(val);
			}
			else
			{
				args->doDxy = 999; // indicates that a string is provided
				args->doDxyStr = strdup(val);
			}
		}

		else if ((strcasecmp("-doNeighborJoining", arv) == 0) || (strcasecmp("-doNJ", arv) == 0))
		{
			args->doNJ = atoi(val);
		}

		else if (strcasecmp("-doDist", arv) == 0)
		{
			args->doDist = atoi(val);
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
		//   --in-jgcd <filename> : joint genotype count distribution input

		else if ((strcasecmp("--in-vcf", arv) == 0) || (strcasecmp("--input", arv) == 0) || (strcasecmp("-i", arv) == 0))
		{
			args->in_vcf_fn = strdup(val);
		}
		else if ((strcasecmp("--in-dm", arv) == 0))
		{
			args->in_dm_fn = strdup(val);
		}
		else if ((strcasecmp("--in-dxy", arv) == 0))
		{
			args->in_dxy_fn = strdup(val);
		}
		else if (strcasecmp("--in-jgcd", arv) == 0)
		{
			args->in_jgcd_fn = strdup(val);
		}

		// TODO if -out has full output name and not prefix, detect it
		// e.g. if you were going to write to a file called "out.amova.csv", and "-out out.amova.csv" is given
		//  then extract prefix from that and use as prefix for other outputs as well
		//  this is a useful feature for snakemake etc
		else if ((strcasecmp("--output", arv) == 0) || (strcasecmp("-out", arv) == 0) || (strcasecmp("-o", arv) == 0))
		{
			args->out_fn = strdup(val);
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
		//   --printJointGenoCountDist/-pJGCD <int> : print joint genotype count distribution
		//   --printJointGenoProbDist/-pJGPD <int> : print joint genotype probability distribution
		//   --printAmovaTable/-pAT <int> : print AMOVA table
		//
		// The following compression levels are available:
		//   0 : no compression
		//   1 : gzip compression
		//   2 : bgzip compression
		//
		//

		// TODO maybe use hypen style here --print-joint-geno-count-dist to be consistent with other commands
		// TODO maybe make <int> optional and choose a default one, most people won't care about this

		else if ((strcasecmp("--printJointGenoCountDist", arv) == 0) || (strcasecmp("--printJGCD", arv) == 0) || (strcasecmp("-pJGCD", arv) == 0))
		{
			args->printJointGenoCountDist = atoi(val);
		}
		else if ((strcasecmp("--printJointGenoProbDist", arv) == 0) || (strcasecmp("--printJGPD", arv) == 0) || (strcasecmp("-pJGPD", arv) == 0))
		{
			args->printJointGenoProbDist = atoi(val);
		}
		else if ((strcasecmp("--printAmovaTable", arv) == 0) || (strcasecmp("--printAT", arv) == 0) || (strcasecmp("-pAT", arv) == 0))
		{
			args->printAmovaTable = atoi(val);
		}
		else if ((strcasecmp("--printDxy", arv) == 0) || (strcasecmp("-pDxy", arv) == 0))
		{
			args->printDxy = atoi(val);
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

		else if (strcasecmp("--region", arv) == 0 || strcmp("-r", arv) == 0)
		{
			args->in_region = strdup(val);
		}

		// N.B. There are two short forms for the regions file command since
		//      -R   easy to remember from bcftools
		//      -rf  people are used to using it in angsd
		else if (strcasecmp("--regions-tab", arv) == 0 || strcasecmp("-rf", arv) == 0 || strcmp("-R", arv) == 0)
		{
			args->in_regions_tab_fn = strdup(val);
		}

		else if (strcasecmp("--regions-bed", arv) == 0 || strcasecmp("-rb", arv) == 0)
		{
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

		else if (strcasecmp("--block-size", arv) == 0 || strcasecmp("-bs", arv) == 0)
		{
			// read block size as float and convert to int
			// to allow for the use of scientific notation (e.g. 1e6)
			args->blockSize = (int) atof(val);
		}

		else if (strcasecmp("--blocks-tab", arv) == 0)
		{
			args->in_blocks_tab_fn = strdup(val);
		}

		else if (strcasecmp("--blocks-bed", arv) == 0)
		{
			args->in_blocks_bed_fn = strdup(val);
		}

		// #####################################
		// #    ARGUMENTS [--long-form/-sf]    #
		// #####################################
		//
		// Use argument commands to specify the parameters to be used in the analyses.
		// Argument commands are of the form `--long-form <int>` or `-sf <int>`
		//   where long form of the argument starts with double dash `--` and separated by hyphen `-`
		//   and short form of the argument starts with single dash `-` and is typically the first letter(s) of the long form
		//

		else if (strcasecmp("--metadata", arv) == 0 || strcasecmp("-m", arv) == 0)
		{
			args->in_mtd_fn = strdup(val);
		}
		else if ((strcasecmp("--verbose", arv) == 0) || (strcasecmp("-v", arv) == 0))
		{
			if (atoi(val) == 0)
			{
				// explicit verbose off, use the default value 0
				fprintf(stderr, "[INFO]\t-> Verbosity disabled explicitly. Will not print any information.\n");
			}
			else if (strIsNumeric(val))
			{
				int verbose_val=atoi(val)-1;
				if(verbose_val>MAX_VERBOSE_LEVEL)
				{
					fprintf(stderr, "\n[WARNING]\tVerbosity level is set to %d, which is greater than the maximum verbosity level %d. Setting verbosity level to %d.\n", verbose_val, MAX_VERBOSE_LEVEL, MAX_VERBOSE_LEVEL);
					verbose_val=MAX_VERBOSE_LEVEL;
				}
				BITSET(VERBOSE, verbose_val);
			}
			else
			{
				NEVER;
			}
			fprintf(stderr, "\n[INFO]\t-> Verbosity level is set to %d.\n", WHICH_BIT_SET1(VERBOSE));
		}

		else if (strcasecmp("--formula", arv) == 0 || strcasecmp("-f", arv) == 0)
		{
			args->formula = strdup(val);
		}

		// [--majorMinorFile/-mmf <filename>]
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

		else if (strcasecmp("--majorMinorFile", arv) == 0 || strcasecmp("-mmf", arv) == 0)
		{
			args->in_majorminor_fn = strdup(val);
		}


		// [--ancDerFile/-adf <filename>]
		//     Filename of the file containing the ancestral and derived alleles for each site in the VCF file.
		//
		//     The file should be in the following format:
		//	   - 1-based indexing for positions [start:included, end:included]
		//	   - Should be sorted by chromosome name and start position.
		//     - Should be tab-delimited.
		//
		// e.g. Chromosome name <TAB> position <TAB> ancestral allele <TAB> derived allele
		//    chr1	100	A	T
		//    chr1	200	A	G
		//
		//    This file can be obtained using ANGSD maf files.
	
		else if (strcasecmp("--ancDerFile", arv) == 0 || strcasecmp("-adf", arv) == 0)
		{
			args->in_ancder_fn = strdup(val);
		}

		//TODO majorminor and ancder file inputs,can we join them?











		else if (strcasecmp("-dev", arv) == 0)
			args->printDev = atoi(val);

		else if (strcasecmp("--isTest", arv) == 0)
			args->isTest = atoi(val);
		else if (strcasecmp("--minInd", arv) == 0)
			args->minInd = atoi(val);

		else if (strcasecmp("--square_distance", arv) == 0)
			args->squareDistance = atoi(val);
		else if (strcasecmp("-sqDist", arv) == 0)
			args->squareDistance = atoi(val);

		else if (strcasecmp("--mThreads", arv) == 0)
			args->mThreads = atoi(val);
		else if (strcasecmp("--mEmIter", arv) == 0)
			args->maxEmIter = atoi(val);

		else if (strcasecmp("--tole", arv) == 0)
			args->tole = atof(val);

		else if (strcasecmp("--gl2gt", arv) == 0)
			args->gl2gt = atoi(val);

		else if (strcasecmp("--seed", arv) == 0)
			args->seed = atoi(val);
		else if (strcasecmp("-seed", arv) == 0)
			args->seed = atoi(val);


		else if ((strcasecmp("-ws", arv) == 0) || (strcasecmp("--windowSize", arv) == 0))
		{
			args->windowSize = (int)atof(val);
		}

		else if ((strcasecmp("-nb", arv) == 0) || (strcasecmp("--nBootstraps", arv) == 0))
		{
			args->nBootstraps = (int)atof(val);
		}
		else if (strcasecmp("-tole", arv) == 0)
			args->tole = atof(val);
		else if (strcasecmp("-isTest", arv) == 0)
			args->isTest = atoi(val);

		else if (strcasecmp("-printMatrix", arv) == 0)
			args->printMatrix = atoi(val);
		else if (strcasecmp("--printMatrix", arv) == 0)
			args->printMatrix = atoi(val);

		else if (strcasecmp("-minInd", arv) == 0)
			args->minInd = atoi(val);
		else if (strcasecmp("-maxIter", arv) == 0)
			args->maxEmIter = atoi(val);
		else if (strcasecmp("-maxEmIter", arv) == 0)
			args->maxEmIter = atoi(val);
		else if (strcasecmp("-mEmIter", arv) == 0)
			args->maxEmIter = atoi(val);
		else if (strcasecmp("-P", arv) == 0)
			args->mThreads = atoi(val);
		else if (strcasecmp("-nThreads", arv) == 0)
			args->mThreads = atoi(val);
		else if (strcasecmp("-gl2gt", arv) == 0)
			args->gl2gt = atoi(val);
		else if (strcasecmp("-h", arv) == 0 || strcasecmp("--help", arv) == 0)
		{
			print_help(stdout);
			exit(0);
		}
		else
		{
			fprintf(stderr, "[ERROR]\tUnknown argument: %s\n", arv);
			exit(1);
		}
		++argv;
	}

	// TODO add 'requires' arg dependency checker, some libs do it like tclap


	if (args->minInd == 0)
	{
		fprintf(stderr, "\n\t-> -minInd 0; will use sites with data for all individuals.\n");
	}
	else if (args->minInd == -1)
	{
		fprintf(stderr, "\n\t-> -minInd not set; will use sites that is nonmissing for both individuals in a pair.\n");
		args->minInd = 2;
	}
	else if (args->minInd == 2)
	{
		fprintf(stderr, "\n\t-> -minInd 2; will use sites that is nonmissing for both individuals in a pair.\n");
	}
	else if (args->minInd == 1 || args->minInd < -1)
	{
		fprintf(stderr, "\n[ERROR]\tMinimum value allowed for minInd is 2.\n");
		exit(1);
	}

	if (args->isTest == 1)
	{
		fprintf(stderr, "Test mode ON\n");
	}

	if (args->out_fn == NULL)
	{
		args->out_fn = strdup("amovaput");
		fprintf(stderr, "\n\t-> -out <output_prefix> not set; will use %s as a prefix for output files.\n", args->out_fn);
	}

	if (args->in_mtd_fn == NULL)
	{
		if (args->doAMOVA != 0)
		{
			fprintf(stderr, "\n[ERROR]\tMust supply -m <Metadata_file> for performing AMOVA.\n");
			exit(1);
		}
	}
	else
	{
		if (args->formula == NULL)
		{
			fprintf(stderr, "\nAMOVA formula not defined, will use all columns in Metadata file %s assuming they are ordered as hierarchical levels.\n", args->in_mtd_fn);
		}
		else
		{
			fprintf(stderr, "\nAMOVA formula is defined as %s; will use the formula to define hierarchical structure in Metadata file %s.\n", args->formula, args->in_mtd_fn);
		}
	}

	if (args->doDist == 1)
	{
		fprintf(stderr, "\n\t-> -doDist is set to 1, will use Dij (1-Sij) dissimilarity index as distance measure.\n");
	}

	if (args->squareDistance == 1)
	{
		fprintf(stderr, "\n\t-> -do_square_distance is set to 1, will square the distance measure (Dij^2).\n");
	}
	else if (args->squareDistance == 0)
	{
		fprintf(stderr, "\n\t-> -do_square_distance is set to 0, will not square the distance measure.\n");
	}
	else
	{
		fprintf(stderr, "\n[ERROR]\t-do_square_distance %d is not available.\n", args->squareDistance);
		exit(1);
	}

	// TODO handle this better
	if (args->in_dm_fn != NULL)
	{
		args->squareDistance = 0;
	}

	if (args->in_vcf_fn == NULL && args->in_dm_fn == NULL)
	{
		fprintf(stderr, "\n[ERROR] Must supply either --in-vcf <VCF_file> or --in-dm <Distance_matrix_file>.\n");
		exit(1);
	}
	else if (args->in_vcf_fn != NULL && args->in_dm_fn != NULL)
	{
		fprintf(stderr, "\n[ERROR] Cannot use --in-vcf %s and --in-dm %s at the same time.\n", args->in_vcf_fn, args->in_dm_fn);
		exit(1);
	}

	if (args->printMatrix != 0)
	{
		fprintf(stderr, "\n[INFO]\t-> -printMatrix %d; will print distance matrix\n", args->printMatrix);
	}

	switch (args->doAMOVA)
	{
	case 0:
	{
		fprintf(stderr, "\n\t-> -doAMOVA is set to 0, will not perform AMOVA.\n");

		if (args->doEM == 1)
		{
			if (args->in_dm_fn != NULL)
			{
				fprintf(stderr, "\n[ERROR] Cannot use -in_dm %s with -doEM 1.\n", args->in_dm_fn);
				exit(1);
			}
			if (args->in_vcf_fn == NULL)
			{
				fprintf(stderr, "\n[ERROR] Must supply -i <input_file> for -doEM 1.\n");
				exit(1);
			}

			fprintf(stderr, "\n\t-> -doEM is set to 1, will use EM algorithm to estimate parameters.\n");
		}

		break;
	}
	case 1: // doAMOVA 1
	{

		if (args->doEM == 0 && args->in_dm_fn == NULL)
		{
			fprintf(stderr, "\n[ERROR]\t-doAMOVA %i requires -doEM 1.\n", args->doAMOVA);
			exit(1);
		}

		if (args->in_dm_fn != NULL)
		{
			fprintf(stderr, "\n-> --in-dm %s is set, will use distance matrix file as data.\n", args->in_dm_fn);
		}
		else
		{

			fprintf(stderr, "\n[INFO]\t-> -doAMOVA 1; will use 10 genotype likelihoods from VCF file GL field.\n");
		}
		break;
	}
	case 2:
	{
		if (args->in_dm_fn != NULL)
		{
			fprintf(stderr, "\n-> --in-dm %s is set, will use distance matrix file as data.\n", args->in_dm_fn);
			args->doAMOVA = 1; // 1: use dm input or gle tag in vcf
		}
		else
		{
			fprintf(stderr, "\n[INFO]\t-> -doAMOVA 2; will use genotypes from VCF file GT field.\n");
		}

		if (args->doEM != 0)
		{
			fprintf(stderr, "\n[ERROR]\t-doAMOVA %i cannot be used with -doEM %i.\n", args->doAMOVA, args->doEM);
			exit(1);
		}
		break;
	}
	case 3:
	{

		if (args->doEM == 0)
		{
			fprintf(stderr, "\n[ERROR]\t-doAMOVA %i requires -doEM 1.\n", args->doAMOVA);
			exit(1);
		}

		if (args->in_dm_fn != NULL)
		{
			fprintf(stderr, "\n-> --in-dm %s is set, will use distance matrix file as data.\n", args->in_dm_fn);
			args->doAMOVA = 1; // 1: use dm input or gle tag in vcf
		}
		else
		{
			fprintf(stderr, "\n[INFO]\t-> -doAMOVA 3; will do both -doAMOVA 1 and -doAMOVA 2.\n");
		}
		break;
	}

	// ---------------------- doAMOVA NOT in {0,1,2,3} ---------------------- //
	default:
	{
		fprintf(stderr, "\n[ERROR]\t-> -doAMOVA %d not recognized\n", args->doAMOVA);
		exit(1);
		break;
	}
	}

	if (args->seed == -1)
	{
		args->seed = time(NULL);
		fprintf(stderr, "\n\t[INFO] -> -seed is not set, will use current time as seed for random number generator: %d.\n", args->seed);
		srand48(args->seed);
	}
	else
	{
		fprintf(stderr, "\n\t-> -seed is set to %d, will use this seed for random number generator.\n", args->seed);
		srand48(args->seed);
	}

	if (args->in_dm_fn != NULL)
	{
		if (args->doEM != 0)
		{
			fprintf(stderr, "\n[ERROR]\t-doEM %i cannot be used with -in_dm %s.\n", args->doEM, args->in_dm_fn);
			exit(1);
		}
	}

	if (args->windowSize == 0)
	{
		fprintf(stderr, "\n[INFO]\t-> -windowSize 0; will not use sliding window\n");
	}
	else
	{
		fprintf(stderr, "\n[INFO]\t-> -windowSize %d; will use sliding windows of size %d\n", args->windowSize, args->windowSize);
	}

	// [dev mode]
#if 1 == DEV
	DEVPRINT("Development mode is on. Will print extra information.\n");
	// BITSET(VERBOSE, 7); // max: 7
#endif

	if (args->doAMOVA > 0)
	{
		IO::requireFile(args->formula, "--formula/-f");
	}

	if (args->doDxy == 0 && args->printDxy == 1)
	{
		fprintf(stderr, "\n[ERROR]\t-> --printDxy 1 requires -doDxy 1.\n");
		exit(1);
	}

	// doNJ 1 requires a distance matrix from -in_dm or -doDist
	// doNJ 2 requires a dxy distance matrix from -in_dxy or -doDxy
	if (args->doNJ == 1 && args->doDist == 0 && args->in_dm_fn == NULL)
	{
		fprintf(stderr, "\n[ERROR]\t-> -doNJ %i requires -doDist 1 or --in-dm <file>.\n", args->doNJ);
		exit(1);
	}
	if (args->doNJ == 2 && args->doDxy == 0 && args->in_dxy_fn == NULL)
	{
		fprintf(stderr, "\n[ERROR]\t-> -doNJ %i requires -doDxy 1 or --in-dxy <file>.\n", args->doNJ);
		exit(1);
	}

	// TODO using regions with block definitions
	// TODO handle empty blocks
	if(args->in_blocks_tab_fn != NULL || args->in_blocks_bed_fn != NULL)
	{
		if(args->in_regions_tab_fn != NULL || args->in_regions_bed_fn != NULL || args->in_region != NULL)
		{
			fprintf(stderr, "\n[ERROR]\tBlock definitions cannot be used with region definitions, yet.\n");
			exit(1);
		}
	}

	return args;
}

void argStruct_destroy(argStruct *args)
{
	FREE(args->in_vcf_fn);
	FREE(args->in_dm_fn);
	FREE(args->in_mtd_fn);
	FREE(args->in_jgcd_fn);
	FREE(args->in_jgpd_fn);
	FREE(args->in_dxy_fn);
	FREE(args->out_fn);
	FREE(args->formula);
	FREE(args->keyCols);
	FREE(args->doDxyStr);

	FREE(args->in_dxy_fn);
	FREE(args->in_region);
	FREE(args->in_regions_tab_fn);
	FREE(args->in_regions_bed_fn);
	FREE(args->in_blocks_tab_fn);
	FREE(args->in_blocks_bed_fn);

	FREE(args->in_majorminor_fn);
	FREE(args->in_ancder_fn);

	delete args;
}

/// @brief argStruct_print - print the arguments to a file pointer
/// @param fp	file pointer
/// @param args	pointer to the argStruct
void argStruct_print(FILE *fp, argStruct *args)
{
	// fprintf(fp, "\nCommand: %s", args->command);//TODO
	fprintf(fp, "\n\t-> -in_vcf_fn %s", args->in_vcf_fn);

	// TODO use lut to store names and values and associatons (e.g. tole maxiter etc assoc with doEM)
	// and if -formula is used, run formulaStruct_get()
	fprintf(fp, "\nngsAMOVA -doAMOVA %d -i %s -out %s -minInd %d -printMatrix %d -m %s -doDist %d -nThreads %d", args->doAMOVA, args->in_vcf_fn, args->out_fn, args->minInd, args->printMatrix, args->in_mtd_fn, args->doDist, args->mThreads);
	if (args->doEM != 0)
	{
		fprintf(fp, " -tole %e ", args->tole);
		fprintf(fp, " -doEM %d ", args->doEM);
		fprintf(fp, " -maxIter %d ", args->maxEmIter);
	}
	if (args->nBootstraps > 0)
	{
		fprintf(fp, " --nBootstraps %d ", args->nBootstraps);
		fprintf(fp, " --blockSize %d ", args->blockSize);
		fprintf(fp, " --seed %d ", args->seed);
	}
	fprintf(fp, "\n");
}
