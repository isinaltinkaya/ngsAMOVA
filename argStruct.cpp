#include "argStruct.h"
#include "paramStruct.h"

argStruct *argStruct_init()
{

	argStruct *args = (argStruct *)calloc(1, sizeof(argStruct));

	args->verbose = 0;

	args->in_vcf_fn = NULL;
	// args->in_sfs_fn = NULL;
	args->in_dm_fn = NULL;
	args->in_mtd_fn = NULL;
	args->out_fn = NULL;

	args->formula = NULL;
	args->keyCols = NULL;
	args->hasColNames = 1;

	args->command = NULL;
	args->blockSize = 0;

	args->windowSize = 0;

	args->doAMOVA = 0;
	args->doEM = 0;

	args->mThreads = 0;
	args->mEmIter = 1e2;

	args->tole = 1e-5;
	args->doTest = 0;

	args->doDist = -1;
	args->do_square_distance = 1;
	

	args->isSim = 0;
	args->isTest = 0;
	args->minInd = -1;

	args->printDev= 0;
	args->printMatrix = 0;
	args->printAmovaTable = 0;
	args->printJointGenoCountDist = 0;
	args->printJointGenoProbDist = 0;


	args->seed = -1;
	args->nBootstraps = 0;

	args->gl2gt = -1;

	return args;
}

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


		if ((strcasecmp("--input", arv) == 0) || (strcasecmp("-in", arv) == 0) || (strcasecmp("-i", arv) == 0)){
			args->in_vcf_fn = strdup(val);
		}
		// else if ((strcasecmp("--inputJointGenoProbDist", arv) == 0) || (strcasecmp("--inputJGPD", arv) == 0) || (strcasecmp("-iJGPD", arv) == 0)){
		// 	args->in_sfs_fn = strdup(val);
		// }
		else if (strcasecmp("-in_dm", arv) == 0)
			args->in_dm_fn = strdup(val);
		else if (strcasecmp("-m", arv) == 0)
			args->in_mtd_fn = strdup(val);
		else if ((strcasecmp("--output", arv) == 0) || (strcasecmp("-out", arv) == 0) || (strcasecmp("-o", arv) == 0)){
			args->out_fn = strdup(val);
		}

		else if (strcasecmp("-dev", arv) == 0)
			args->printDev = atoi(val);


		else if ((strcasecmp("--printJointGenoCountDist", arv) == 0) || (strcasecmp("--printJGCD", arv) == 0) || (strcasecmp("-pJGCD", arv) == 0)){
			args->printJointGenoCountDist = atoi(val);
		}
		else if ((strcasecmp("--printJointGenoProbDist", arv) == 0) || (strcasecmp("--printJGPD", arv) == 0) || (strcasecmp("-pJGPD", arv) == 0)){
			args->printJointGenoProbDist = atoi(val);
		}
		else if ((strcasecmp("--printAmovaTable", arv) == 0) || (strcasecmp("--printAT", arv) == 0) || (strcasecmp("-pAT", arv) == 0)){
			args->printAmovaTable = atoi(val);
		}

		else if ((strcasecmp("--verbose", arv) == 0) || (strcasecmp("-v", arv) == 0)) 
		{
			if (val == NULL)
				args->verbose = 1;
			else if (isdigit(val[0]))
				args->verbose = atoi(val);
			else
				args->verbose = 1;
		}

		else if (strcasecmp("--isSim", arv) == 0)
			args->isSim = atoi(val);
		else if (strcasecmp("--isTest", arv) == 0)
			args->isTest = atoi(val);
		else if (strcasecmp("--minInd", arv) == 0)
			args->minInd = atoi(val);

		else if (strcasecmp("--doDist", arv) == 0)
			args->doDist = atoi(val);
		else if (strcasecmp("--do_square_distance", arv) == 0)
			args->do_square_distance = atoi(val);

		else if (strcasecmp("--doAMOVA", arv) == 0)
			args->doAMOVA = atoi(val);
		else if (strcasecmp("--doEM", arv) == 0)
			args->doEM = atoi(val);

		else if (strcasecmp("--mThreads", arv) == 0)
			args->mThreads = atoi(val);
		else if (strcasecmp("--mEmIter", arv) == 0)
			args->mEmIter = atoi(val);

		else if (strcasecmp("--tole", arv) == 0)
			args->tole = atof(val);
		else if (strcasecmp("--doTest", arv) == 0)
			args->doTest = atoi(val);

		else if (strcasecmp("--gl2gt", arv) == 0)
			args->gl2gt = atoi(val);

		else if (strcasecmp("--seed", arv) == 0)
			args->seed = atoi(val);


		// read block size as float and convert to int
		// this is to allow for the use of scientific notation (e.g. 1e6)
		else if (strcasecmp("-bs", arv) == 0)
			args->blockSize = (int) atof(val);
		else if (strcasecmp("--blockSize", arv) == 0)
			args->blockSize = (int) atof(val);

		else if (strcasecmp("-ws", arv) == 0)
			args->windowSize = (int) atof(val);
		else if (strcasecmp("--windowSize", arv) == 0)
			args->windowSize = (int) atof(val);

		else if (strcasecmp("-bSize", arv) == 0)
			args->blockSize = (int) atof(val);
		else if (strcasecmp("-nb", arv) == 0)
			args->nBootstraps = (int) atof(val);


		else if (strcasecmp("-f", arv) == 0)
		{
			args->formula = strdup(val);
		}
		else if (strcasecmp("--formula", arv) == 0)
		{
			args->formula = strdup(val);
		}
		else if (strcasecmp("--hasColNames", arv) == 0)
			args->hasColNames = atoi(val);
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
		else if (strcasecmp("-doDist", arv) == 0)
			args->doDist = atoi(val);
		else if (strcasecmp("-sqDist", arv) == 0)
			args->do_square_distance = atoi(val);
		else if (strcasecmp("-minInd", arv) == 0)
			args->minInd = atoi(val);
		else if (strcasecmp("-doTest", arv) == 0)
			args->doTest = atoi(val);
		else if (strcasecmp("-maxIter", arv) == 0)
			args->mEmIter = atoi(val);
		else if (strcasecmp("-maxEmIter", arv) == 0)
			args->mEmIter = atoi(val);
		else if (strcasecmp("-mEmIter", arv) == 0)
			args->mEmIter = atoi(val);
		else if (strcasecmp("-P", arv) == 0)
			args->mThreads = atoi(val);
		else if (strcasecmp("-nThreads", arv) == 0)
			args->mThreads = atoi(val);
		else if (strcasecmp("-gl2gt", arv) == 0)
			args->gl2gt = atoi(val);
		else if (strcasecmp("-h", arv) == 0 || strcasecmp("--help", arv) == 0)
		{
			free(args);
			usage(stdout);
			exit(0);
		}
		else
		{
			fprintf(stderr, "[ERROR]\tUnknown argument: %s\n", arv);
			exit(1);
		}
		++argv;
	}
	// TODO add 'requires' arg dependency checker

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

	if (args->doAMOVA != 3 && args->doTest == 1)
	{
		fprintf(stderr, "\n[ERROR]\t-doTest 1 requires -doAMOVA 3.\n");
		exit(1);
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
			if (args->hasColNames == 0)
			{
				fprintf(stderr, "\n[ERROR]\tAMOVA formula is defined but -hasColnames is set to 0. Metadata file must have column names if formula is defined.\n");
			}
		}
	}

	if (args->doDist == 1)
	{
		fprintf(stderr, "\n\t-> -doDist is set to 1, will use Dij (1-Sij) dissimilarity index as distance measure.\n");
	}
	else
	{
		fprintf(stderr, "\n[ERROR]\t-doDist %d is not available.\n", args->doDist);
		exit(1);
	}

	if (args->do_square_distance == 1)
	{
		fprintf(stderr, "\n\t-> -do_square_distance is set to 1, will square the distance measure (Dij^2).\n");
	}
	else if (args->do_square_distance == 0)
	{
		fprintf(stderr, "\n\t-> -do_square_distance is set to 0, will not square the distance measure.\n");
	}else{
		fprintf(stderr, "\n[ERROR]\t-do_square_distance %d is not available.\n", args->do_square_distance);
		exit(1);
	}

	if (args->in_vcf_fn == NULL && args->in_dm_fn == NULL)
	{
		fprintf(stderr, "\n[ERROR] Must supply either -in <VCF_file> or -in_dm <Distance_matrix_file>.\n");
		exit(1);
	}else if (args->in_vcf_fn != NULL && args->in_dm_fn != NULL)
	{
		fprintf(stderr, "\n[ERROR] Cannot use -in %s with -in_dm %s.\n", args->in_vcf_fn, args->in_dm_fn);
		exit(1);
	}

	if (args->printMatrix == 1){
		fprintf(stderr, "\n[INFO]\t-> -printMatrix 1; will print distance matrix\n");
	}else if (args->printMatrix == 0){
		// fprintf(stderr, "\n[INFO]\t-> -printMatrix 0; will not print distance matrix\n");
	}else{
		fprintf(stderr, "\n[ERROR]\t-> -printMatrix %d is not available.\n", args->printMatrix);
		exit(1);
	}

	switch (args->doAMOVA)
	{
		case 0:
		{
			fprintf(stderr, "\n\t-> -doAMOVA is set to 0, will not perform AMOVA.\n");

			if (args->doEM == 0)
			{
				fprintf(stderr, "\n\t-> Nothing to do.\n");
				exit(1);
			}
			else if (args->doEM == 1)
			{
				if (args->in_dm_fn != NULL)
				{
					fprintf(stderr, "\n[ERROR] Cannot use -in_dm %s with -doEM 1.\n", args->in_dm_fn);
					exit(1);
				}
				if (args->in_vcf_fn == NULL)
				{
					fprintf(stderr, "\n[ERROR] Must supply -in <input_file> for -doEM 1.\n");
					exit(1);
				}

				fprintf(stderr, "\n\t-> -doEM is set to 1, will use EM algorithm to estimate parameters.\n");
			}

			break;
		}
		case 1:
		{

			if (args->doEM == 0 && args->in_dm_fn == NULL)
			{
				fprintf(stderr, "\n[ERROR]\t-doAMOVA %i requires -doEM 1.\n", args->doAMOVA);
				exit(1);
			}

			if (args->in_dm_fn != NULL)
			{
				fprintf(stderr, "\n-> -in_dm %s is set, will use distance matrix file as data.\n", args->in_dm_fn);
			}
			else
			{

				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 1; will use 10 genotype likelihoods from VCF GL tag.\n");
			}
			break;
		}
		case 2:
		{
			if (args->in_dm_fn != NULL)
			{
				fprintf(stderr, "\n-> -in_dm %s is set, will use distance matrix file as data.\n", args->in_dm_fn);
				args->doAMOVA = 1; // 1: use dm input or gle tag in vcf
			}
			else
			{

				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 2; will use genotypes from VCF GT tag.\n");
			}
			
			if (args->doEM != 0){
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
				fprintf(stderr, "\n-> -in_dm %s is set, will use distance matrix file as data.\n", args->in_dm_fn);
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


	if(args->nBootstraps>0){

		fprintf(stderr, "\n\t-> -nBootstraps %d is set, will perform %d bootstraps for AMOVA significance testing.\n", args->nBootstraps, args->nBootstraps);
		if(args->blockSize==0){
			fprintf(stderr, "\n[ERROR] -blockSize must be set to a positive integer when -nBootstraps is set.\n");
			exit(1);
		}else{
			fprintf(stderr, "\n\t-> -blockSize %d is set, will use %d as the genomic block size for block bootstrapping.\n", args->blockSize, args->blockSize);

			if(args->seed == -1){
				args->seed = time(NULL);
				fprintf(stderr, "\n\t[INFO] -> -seed is not set, will use current time as seed for random number generator: %d.\n", args->seed);
				srand48(args->seed);
			}else{
				fprintf(stderr, "\n\t-> -seed is set to %d, will use this seed for random number generator.\n", args->seed);
				srand48(args->seed);
			}
		}

	}else if(args->nBootstraps<0){
		fprintf(stderr, "\n[ERROR]\t-> -nBootstraps should be a positive integer or 0. You entered a negative value: %d.\n", args->nBootstraps);
		exit(1);
	}else{
		if(args->blockSize!=0){
			fprintf(stderr, "\n[ERROR] -blockSize is set to %d, but -nBootstraps is not set. Define both to perform block bootstrapping.\n", args->blockSize);
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


	if(args->windowSize == 0){
		fprintf(stderr, "\n[INFO]\t-> -windowSize 0; will not use sliding window\n");
	}else{
		fprintf(stderr, "\n[INFO]\t-> -windowSize %d; will use sliding windows of size %d\n", args->windowSize, args->windowSize);
	}


	return args;
}

void argStruct_destroy(argStruct *args)
{
	FREE(args->in_vcf_fn);
	FREE(args->in_jgcd_fn);
	FREE(args->in_jgpd_fn);
	FREE(args->in_dm_fn);
	FREE(args->in_mtd_fn);
	FREE(args->out_fn);
	FREE(args->formula);
	FREE(args->keyCols);
	FREE(args);
}

/// @brief argStruct_print - print the arguments to a file pointer
/// @param fp	file pointer
/// @param args	pointer to the argStruct
void argStruct_print(FILE* fp, argStruct *args){
	// fprintf(fp, "\nCommand: %s", args->command);//TODO
	fprintf(fp, "\n\t-> -in_vcf_fn %s", args->in_vcf_fn);

	// TODO print based on analysis type, and to args file, collect all from args automatically
	fprintf(fp, "\nngsAMOVA -doAMOVA %d -doTest %d -in %s -out %s -isSim %d -minInd %d -printMatrix %d -m %s -doDist %d -maxIter %d -nThreads %d", args->doAMOVA, args->doTest, args->in_vcf_fn, args->out_fn, args->isSim, args->minInd, args->printMatrix, args->in_mtd_fn, args->doDist, args->mEmIter, args->mThreads);
	if (args->doEM != 0)
	{
		fprintf(fp, " -tole %e ", args->tole);
		fprintf(fp, " -maxIter %d ", args->mEmIter);
	}
	fprintf(fp, "\n");
}

/// @brief check_consistency_args_pars - check consistency between arguments and parameters
/// @param args pointer to argStruct
/// @param pars pointer to paramStruct
void check_consistency_args_pars(argStruct *args, paramStruct *pars){

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
		
			if(pars->in_ft == IN_DM && args->printMatrix == 1){
				fprintf(stderr, "\n\n[ERROR]\tCannot print distance matrix since input file is already a distance matrix; will exit!\n\n");
				exit(1);
			}

}