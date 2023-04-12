#include "argStruct.h"
#include "paramStruct.h"
#include "io.h"

// default: 0 (verbose mode off)
u_char VERBOSE = 0;

argStruct *argStruct_init()
{

	argStruct *args = (argStruct *)malloc(sizeof(argStruct));

	args->in_vcf_fn = NULL;
	args->in_dm_fn = NULL;
	args->in_mtd_fn = NULL;
	args->in_jgcd_fn = NULL;
	args->in_jgpd_fn = NULL;
	args->in_blb_fn = NULL;
	args->in_dxy_fn = NULL;

	args->out_fn = NULL;

	args->formula = NULL;
	args->keyCols = NULL;

	args->command = NULL;
	args->blockSize = 0;

	args->windowSize = 0;

	args->doAMOVA = 0;
	args->doEM = 0;
	args->doDxy = 0;
	args->doDxyStr = NULL;
	args->doNJ = 0;

	args->mThreads = 0;
	args->maxEmIter = 100;

	args->tole = 1e-5;

	args->doDist = -1;
	args->square_distance = 1;

	args->isSim = 0;
	args->isTest = 0;
	args->minInd = -1;

	args->printDev = 0;
	args->printMatrix = 0;
	args->printAmovaTable = 0;
	args->printJointGenoCountDist = 0;
	args->printJointGenoProbDist = 0;
	args->printDxy = 0;

	args->seed = -1;
	args->nBootstraps = 0;

	args->gl2gt = -1;

	return args;
}

// TODO check multiple of same argument
/// @brief argStruct_get read command line arguments
/// @param argc
/// @param argv
/// @return pointer to argStruct
argStruct *argStruct_get(int argc, char **argv)
{

	argStruct *args = argStruct_init();

	while (*argv)
	{

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
		// Action commands are of the form `-doXXXX <value>`
		//   where XXXX defines the general type of analysis to be performed
		//   and <value> specifies the exact analysis to be performed.
		// For example, `-doAMOVA <value>` specifies to perform AMOVA analysis.
		//   and <value> specifies the exact AMOVA analysis to be performed.
		//   e.g. `-doAMOVA 1` specifies to perform AMOVA analysis with the genotype likelihoods.
		//
		// The following action commands are available:
		//   -doAMOVA <value> : perform AMOVA analysis
		//   -doEM <value> : perform EM optimization
		//   -doDxy <value> : estimate Dxy
		//   -doNJ <value> : do neighbor-joining
		//   -doDist <value> : calculate pairwise distances

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

		// TODO this is prefix, but add output name checker to make it play nicer with snakemake
		// e.g. --output out.amova.csv should detect that user meant to use out as prefix not the whole thing
		// if one of the expected output filename is out.amova.csv
		else if ((strcasecmp("--output", arv) == 0) || (strcasecmp("-out", arv) == 0) || (strcasecmp("-o", arv) == 0))
		{
			args->out_fn = strdup(val);
		}

		// #################################################################
		// #    PRINTING COMMANDS [--printXxXxx/--printXX/-pXX <value>]    #
		// #################################################################
		//
		// Use printing commands to specify the output files to be generated.
		// This is only needed for output files that are not the default output files
		//   associated with the analyses specified by the action commands.
		//
		// Printing commands are of the form `--printXxXxx <value>` or `--printXX <value>` or `-pXX <value>`
		//   where XxXxx is the long form of the file type to be printed,
		//   XX is the short form of the file type (typically  the first letter(s) of the long form),
		//   and <value> defines the compression level of the output file.
		// Output file names are automatically generated based on the value of the `--output` argument
		//   and the file type to be printed.
		//
		// The following printing commands are available:
		//   --printJointGenoCountDist/-pJGCD <value> : print joint genotype count distribution
		//   --printJointGenoProbDist/-pJGPD <value> : print joint genotype probability distribution
		//   --printAmovaTable/-pAT <value> : print AMOVA table
		//
		// The following compression levels are available:
		//   0 : no compression
		//   1 : gzip compression
		//   2 : bgzip compression
		//
		//

		// TODO maybe use hypen style here --print-joint-geno-count-dist to be consistent with other commands
		// TODO maybe make <value> optional and choose a default one, most people won't care about this

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

		// #####################################
		// #    ARGUMENTS [--long-form/-sf]    #
		// #####################################
		//
		// Use argument commands to specify the parameters to be used in the analyses.
		// Argument commands are of the form `--long-form <value>` or `-sf <value>`
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
				BITSET(VERBOSE, (atoi(val) - 1));
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

		else if ((strcasecmp("--block_bed", arv) == 0) || (strcasecmp("-bf", arv) == 0))
		{
			args->in_blb_fn = strdup(val);
		}

		else if (strcasecmp("-dev", arv) == 0)
			args->printDev = atoi(val);

		else if (strcasecmp("--isSim", arv) == 0)
			args->isSim = atoi(val);
		else if (strcasecmp("--isTest", arv) == 0)
			args->isTest = atoi(val);
		else if (strcasecmp("--minInd", arv) == 0)
			args->minInd = atoi(val);

		else if (strcasecmp("--square_distance", arv) == 0)
			args->square_distance = atoi(val);
		else if (strcasecmp("-sqDist", arv) == 0)
			args->square_distance = atoi(val);

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

		// read block size as float and convert to int
		// this is to allow for the use of scientific notation (e.g. 1e6)
		else if ((strcasecmp("-bs", arv) == 0) || (strcasecmp("--blockSize", arv) == 0))
		{
			args->blockSize = (int)atof(val);
		}

		else if ((strcasecmp("-ws", arv) == 0) || (strcasecmp("--windowSize", arv) == 0))
		{
			args->windowSize = (int)atof(val);
		}

		else if ((strcasecmp("-nb", arv) == 0) || (strcasecmp("--nBootstraps", arv) == 0))
		{
			args->nBootstraps = (int)atof(val);
		}
		else if (strcasecmp("-seed", arv) == 0)
			args->seed = atoi(val);
		else if (strcasecmp("-doAMOVA", arv) == 0)
			args->doAMOVA = atoi(val);
		else if (strcasecmp("-doEM", arv) == 0)
			args->doEM = atoi(val);
		else if (strcasecmp("-tole", arv) == 0)
			args->tole = atof(val);
		else if (strcasecmp("-isSim", arv) == 0)
			args->isSim = atoi(val);
		else if (strcasecmp("-isTest", arv) == 0)
			args->isTest = atoi(val);

		else if (strcasecmp("-printMatrix", arv) == 0)
			args->printMatrix = atoi(val);
		else if (strcasecmp("--printMatrix", arv) == 0)
			args->printMatrix = atoi(val);

		else if (strcasecmp("-doDist", arv) == 0)
			args->doDist = atoi(val);
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
			free(args);
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

	if (args->isSim > 1 || args->isSim < 0)
	{
		fprintf(stderr, "\n[ERROR]\tArgument isSim is set to %d\n", args->isSim);
		exit(1);
	}

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
	else
	{
		if (args->in_dm_fn == NULL)
		{
			fprintf(stderr, "\n[ERROR]\t-doDist %d is not available.\n", args->doDist);
			exit(1);
		}
	}

	if (args->square_distance == 1)
	{
		fprintf(stderr, "\n\t-> -do_square_distance is set to 1, will square the distance measure (Dij^2).\n");
	}
	else if (args->square_distance == 0)
	{
		fprintf(stderr, "\n\t-> -do_square_distance is set to 0, will not square the distance measure.\n");
	}
	else
	{
		fprintf(stderr, "\n[ERROR]\t-do_square_distance %d is not available.\n", args->square_distance);
		exit(1);
	}

	// TODO handle this better
	if (args->in_dm_fn != NULL)
	{
		args->square_distance = 0;
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

	if (args->nBootstraps > 0)
	{

		fprintf(stderr, "\n\t-> --nBootstraps %d is set, will perform %d bootstraps for AMOVA significance testing.\n", args->nBootstraps, args->nBootstraps);
		if (args->blockSize == 0)
		{
			fprintf(stderr, "\n[ERROR] -blockSize must be set to a positive integer when --nBootstraps is set.\n");
			exit(1);
		}
		else
		{
			fprintf(stderr, "\n\t-> -blockSize %d is set, will use %d as the genomic block size for block bootstrapping.\n", args->blockSize, args->blockSize);

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
		}
	}
	else if (args->nBootstraps < 0)
	{
		fprintf(stderr, "\n[ERROR]\t-> --nBootstraps should be a positive integer or 0. You entered a negative value: %d.\n", args->nBootstraps);
		exit(1);
	}
	else
	{
		if (args->blockSize != 0)
		{
			fprintf(stderr, "\n[ERROR] --blockSize is set to %d, but --nBootstraps is not set. Define both to perform block bootstrapping.\n", args->blockSize);
			exit(1);
		}
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

	return args;
}

void argStruct_destroy(argStruct *args)
{
	FREE(args->in_vcf_fn);
	FREE(args->in_dm_fn);
	FREE(args->in_mtd_fn);
	FREE(args->in_jgcd_fn);
	FREE(args->in_jgpd_fn);
	FREE(args->in_blb_fn);
	FREE(args->in_dxy_fn);
	FREE(args->out_fn);
	FREE(args->formula);
	FREE(args->keyCols);
	FREE(args->doDxyStr);
	FREE(args);
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
	fprintf(fp, "\nngsAMOVA -doAMOVA %d -i %s -out %s -isSim %d -minInd %d -printMatrix %d -m %s -doDist %d -nThreads %d", args->doAMOVA, args->in_vcf_fn, args->out_fn, args->isSim, args->minInd, args->printMatrix, args->in_mtd_fn, args->doDist, args->mThreads);
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

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(argStruct *args, paramStruct *pars)
{

	if (args->minInd == pars->nInd)
	{
		fprintf(stderr, "\n\t-> -minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 (all).\n", args->minInd, pars->nInd);
		args->minInd = 0;
	}

	if (pars->nInd == 1)
	{
		fprintf(stderr, "\n\n[ERROR]\tOnly one sample; will exit\n\n");
		exit(1);
	}

	if (pars->nInd < args->minInd)
	{
		fprintf(stderr, "\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n", args->minInd, pars->nInd);
		exit(1);
	}

	if (pars->in_ft == IN_DM && args->printMatrix == 1)
	{
		fprintf(stderr, "\n\n[ERROR]\tCannot print distance matrix since input file is already a distance matrix; will exit!\n\n");
		exit(1);
	}
}