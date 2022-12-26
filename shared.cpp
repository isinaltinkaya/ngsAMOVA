/*
 *
 * [Parameters]
 *
 * idea from angsd funkyPars
 * at https://github.com/ANGSD/angsd/blob/master/shared.cpp
 *
 */

#include "shared.h"
#include <stdlib.h>

using size_t = decltype(sizeof(int));

//
// 0 1 2
// 00 01 02
// MMMM MMMm MMmm
//
// 3 4 5
// 10 11 12
// MmMM MmMm Mmmm
//
// 6 7 8
// 20 21 22
// mmMM mmMm mmmm
extern const int get_3x3_idx[3][3] = {
	{0, 1, 2},
	{3, 4, 5},
	{6, 7, 8}};

char *get_time()
{
	time_t current_time;
	struct tm *local_time;
	current_time = time(NULL);
	local_time = localtime(&current_time);
	return (asctime(local_time));
}

FILE *IO::getFILE(const char *fname, const char *mode)
{
	if (strcmp(mode, "r") == 0)
	{
		fprintf(stderr, "\n\t-> Reading file: %s\n", fname);
	}
	FILE *fp;
	if (NULL == (fp = fopen(fname, mode)))
	{
		fprintf(stderr, "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n", __FILE__, __FUNCTION__, fname);
		exit(1);
	}
	return fp;
}

char *IO::setFileName(const char *a, const char *b)
{
	char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
	strcpy(c, a);
	strcat(c, b);
	// fprintf(stderr,"\t-> Opening output file for writing: %s\n",c);
	return c;
}

FILE *IO::openFILE(char *c)
{
	fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
	FILE *fp = getFILE(c, "w");
	return fp;
}

FILE *IO::openFILE(const char *a, const char *b)
{
	char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
	strcpy(c, a);
	strcat(c, b);
	fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
	FILE *fp = getFILE(c, "w");
	free(c);
	return fp;
}

int IO::inspectFILE::count_nColumns(char *line, const char *delims)
{

	char *str = NULL;
	str = strdup(line);

	char *p = NULL;
	int i = 0;
	p = strtok(str, delims);
	while (p != NULL)
	{
		i++;
		p = strtok(NULL, delims);
	}

	free(p);
	free(str);
	return i;
}

int IO::readFILE::SFS(FILE *in_sfs_ff, const char *delims, DATA::samplesStruct *SAMPLES)
{

	char sfs_buf[1024];
	while (fgets(sfs_buf, 1024, in_sfs_ff))
	{

		char *tok = strtok(sfs_buf, delims);
		char *col = tok;

		for (int coli = 0; coli < 9 - 1; coli++)
		{
			tok = strtok(NULL, "\t \n");
			col = tok;
		}

#if 1
		fprintf(stderr, "\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "--------------------------------------------------");
		fprintf(stderr, "\n");

		fprintf(stderr, "%s", col);
		fprintf(stderr, "\n");
		fprintf(stderr, "--------------------------------------------------");
		fprintf(stderr, "\n");
#endif
	}

	return 0;
}

int IO::readFILE::METADATA(DATA::metadataStruct *MTD, FILE *in_mtd_ff, int whichCol, const char *delims, DATA::samplesStruct *SAMPLES)
{

	fprintf(stderr, "\n");
	fprintf(stderr, "--------------------------------------------------");
	fprintf(stderr, "\n");

	// int whichCol=args->whichCol;
	int nCols = 0;
	int checkCol = 1;

	// TODO map is probably better due to sorting issues. we cannot have all strata sorted
	char mt_buf[1024];

	while (fgets(mt_buf, 1024, in_mtd_ff))
	{

		if (checkCol == 1)
		{
			nCols = IO::inspectFILE::count_nColumns(mt_buf, delims);
			fprintf(stderr, "\n\t-> Number of columns in input metadata file: %d\n", nCols);
			checkCol = 0;

			if (whichCol != -1 && whichCol > nCols)
			{
				fprintf(stderr, "\n[ERROR]\tColumn %d was chosen, but input metadata file only contains %d columns; will exit!\n\n", whichCol, nCols);
				exit(1);
			}
		}

		// TODO strtok_r? do we need thread safety here?

		char *tok = strtok(mt_buf, delims);
		char *group_id = tok;

		for (int coli = 0; coli < whichCol - 1; coli++)
		{
			tok = strtok(NULL, "\t \n");
			group_id = tok;
		}

		// increase the size of Strata
		if (MTD->nStrata > (int)MTD->_strataArr)
		{
			fprintf(stderr, "->->->increase the size of Strata S[]!! Found %d MTD->nStrata and %d MTD->_strataArr\n", (int)MTD->nStrata, (int) MTD->_strataArr);
		}

		// if not the first loop
		if (MTD->strataArr[MTD->nStrata].id != NULL)
		{

			// group id changed
			if (strcmp(MTD->strataArr[MTD->nStrata].id, group_id) != 0)
			{
				MTD->nStrata++;
				MTD->strataArr[MTD->nStrata].id = strdup(group_id);
			}

			SAMPLES->sampleArr[MTD->nInds_total] += (int)pow(2, (int)MTD->nStrata);

			MTD->strataArr[MTD->nStrata].nInds++;
			MTD->nInds_total++;
		}
		else
		{
			// if first loop

			// set bit
			// nth bit is associated with strata_n
			SAMPLES->sampleArr[MTD->nInds_total] += (int)pow(2, (int)MTD->nStrata);

			MTD->strataArr[MTD->nStrata].id = strdup(group_id);
			MTD->strataArr[MTD->nStrata].nInds++;
			MTD->nInds_total++;
		}

		// TODO then plug in all pairs associated with ind if ind==indid in header in loop
	}
	// add the last one since increase is at the beginning of change
	MTD->nStrata++;

#if 0
	for (int sti=0; sti<MTD->nStrata; sti++){
		fprintf(stderr,"\n-> Strata %s contains %d individuals.",MTD->strataArr[sti].id,MTD->strataArr[sti].nInds);
		fprintf(stderr,"\n-> Individual indexes are:\t");
		for(int ii=0; ii<MTD->nInds_total;ii++){
			if( (SAMPLES->sampleArr[ii] & (1 << sti))){
				fprintf(stderr,"%i",ii);
			}
		}
		fprintf(stderr,"\n");
	}

	fprintf(stderr,"\n");
	fprintf(stderr,"--------------------------------------------------");
	fprintf(stderr,"\n");
#endif

	return 0;
}

/* IO::print::Array
 * Prints the elements of an array to a file or stream.
 *
 * @param arr: the array to print
 * @param N: the number of rows, dimension 1
 * @param M: the number of columns, dimension 2
 * @param out: the file or stream to print to
 * @param sep: the character to use as a separator
 *
 * @example
 * IO::print::Array(stdout, myArray, DOUBLE, 3, 3, ',');
 */
void IO::print::Array(FILE *fp, double *arr, size_t N, size_t M, char sep)
{

	ASSERT(arr != NULL);

	double *p = arr;
	size_t n = 0;
	size_t m = 0;

	while (n < N && p < arr + N * M)
	{

		fprintf(fp, "%f%c", *p, (n == M - 1 && m == M - 1) ? '\n' : sep);

		m++;

		// If the column index m reaches the end of the row, reset it and increment the row index n
		if (m == M)
		{
			m = 0;
			n++;
		}

		// Increment the pointer to the next element of the array
		p++;
	}
}

// IO::print::Array
// :overload: int array
void IO::print::Array(FILE *fp, int *arr, size_t N, size_t M, char sep)
{

	ASSERT(arr != NULL);

	int *p = arr;
	size_t n = 0;
	size_t m = 0;

	while (n < N && p < arr + N * M)
	{

		fprintf(fp, "%d%c", *p, (n == M - 1 && m == M - 1) ? '\n' : sep);

		m++;

		if (m == M)
		{
			m = 0;
			n++;
		}

		p++;
	}
}

// Check if file exists
// @param in_fn	input filename
// @return		1 if file exists; 0 otherwise
// credit: angsd/aio.cpp
int file_exists(const char *in_fn)
{
	struct stat buffer;
	return (stat(in_fn, &buffer) == 0);
}

void usage(FILE *fp)
{
	// fprintf(stderr,"");
	// fprintf(stderr,"\n");
	// fprintf(stderr,"  --help         : Print this help\n");
	// fprintf(stderr,"\t--in\t\t\t: input VCF/BCF filed\n");

	fprintf(fp,
			"\n"
			"Program: ngsAMOVA\n");
	// "Version: %s (using htslib %s)\n\n", program_version(), hts_version());
	// "build(%s %s)\n",__DATE__,__TIME__);

	fprintf(fp,
			"\n"
			"Usage:\tngsAMOVA <command> [options]\n"
			"\n"
			"Tool for performing the Analysis of Molecular Variance [AMOVA]\n"
			"\n"
			"Commands:\n"
			"\t-- Analyses\n"
			"\n"
			"Options:\n"
			" -in/-i\n"
			"\n"
			"\t-s\n"
			"\t-m\n"
			"\t-out/o\n"
			"\n"

			"\t-bs/bSize\n"
			"\t-mCol\n"
			"\t-seed\n"
			"\t-doAMOVA\n"
			"\t-doInd\n"
			"\t-ind1\n"
			"\t-ind2\n"
			"\t-printMatrix\n"
			"\t-doDist\n"
			"\t-sqDist (default: 1)\n"
			"\t-minInd\n"
			"\t-doTest\n"
			"\t-maxIter/maxEmIter/mEmIter\n"
			"\t-P/nThreads (default: 1)\n"
			"\t-gl2gt\n"
			"\n");

	//
	//
	// exit(0);
}

argStruct *argStruct_init()
{

	argStruct *args = (argStruct *)calloc(1, sizeof(argStruct));

	args->in_fn = NULL;
	args->in_sfs_fn = NULL;
	args->in_mtd_fn = NULL;
	args->out_fp = NULL;

	args->whichCol = -1;

	args->blockSize = 0;

	args->seed = -1;
	args->doAMOVA = 0;

	args->mThreads = 0;

	args->mEmIter = 1e2;

	args->tole = 1e-10;

	args->doTest = 0;

	args->doDist = -1;
	args->sqDist = 1;

	args->isSim = 0;
	args->minInd = -1;

	args->printMatrix = 0;

	args->doInd = 0;
	args->ind1 = -1;
	args->ind2 = -1;

	args->gl2gt = -1;

	return args;
}

//
// void *argStruct_destroy(argStruct *args){
// if(args->in_fn){
// free(args->in_fn);
// args->in_fn=NULL;
// }
//
// }
//

argStruct *argStruct_get(int argc, char **argv)
{

	argStruct *args = argStruct_init();

	while (*argv)
	{

		char *arv = *argv;
		char *val = *(++argv);

		if (strcasecmp("-in", arv) == 0)
			args->in_fn = strdup(val);
		else if (strcasecmp("-i", arv) == 0)
			args->in_fn = strdup(val);
		else if (strcasecmp("-s", arv) == 0)
			args->in_sfs_fn = strdup(val);
		else if (strcasecmp("-m", arv) == 0)
			args->in_mtd_fn = strdup(val);
		else if (strcasecmp("-out", arv) == 0)
			args->out_fp = strdup(val);
		else if (strcasecmp("-o", arv) == 0)
			args->out_fp = strdup(val);
		else if (strcasecmp("-bs", arv) == 0)
			args->blockSize = atoi(val);
		else if (strcasecmp("-bSize", arv) == 0)
			args->blockSize = atoi(val);
		else if (strcasecmp("-mCol", arv) == 0)
			args->whichCol = atoi(val);
		else if (strcasecmp("-seed", arv) == 0)
			args->seed = atoi(val);
		else if (strcasecmp("-doAMOVA", arv) == 0)
			args->doAMOVA = atoi(val);
		else if (strcasecmp("-doInd", arv) == 0)
			args->doInd = atoi(val);
		else if (strcasecmp("-ind1", arv) == 0)
			args->ind1 = atoi(val);
		else if (strcasecmp("-ind2", arv) == 0)
			args->ind2 = atoi(val);
		else if (strcasecmp("-tole", arv) == 0)
			args->tole = atof(val);
		else if (strcasecmp("-isSim", arv) == 0)
			args->isSim = atoi(val);
		else if (strcasecmp("-printMatrix", arv) == 0)
			args->printMatrix = atoi(val);
		else if (strcasecmp("-doDist", arv) == 0)
			args->doDist = atoi(val);
		else if (strcasecmp("-sqDist", arv) == 0)
			args->sqDist = atoi(val);
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
			return 0;
		}
		else
		{
			fprintf(stderr, "Unknown arg:%s\n", arv);
			free(args);
			return 0;
		}
		++argv;
	}

	if (args->isSim > 1 || args->isSim < 0)
	{
		fprintf(stderr, "\n[ERROR]\tArgument isSim is set to %d\n", args->isSim);
		free(args);
		return 0;
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
	else if (args->minInd == 1)
	{
		fprintf(stderr, "\n[ERROR]\tMinimum value allowed for minInd is 2; will exit!\n");
		free(args);
		return 0;
	}

	if (args->in_fn == NULL)
	{
		fprintf(stderr, "\n[ERROR]\tMust supply -in <input_file>; will exit!\n");
		free(args);
		return 0;
	}

	// if (args->isSim == 0)
	// {
	// 	fprintf(stderr, "Must use -isSim 1 for now\n");
	// 	free(args);
	// 	return 0;
	// }

	if (args->out_fp == NULL)
	{
		args->out_fp = strdup("amovaput");
		fprintf(stderr, "\n\t-> -out <output_prefix> not set; will use %s as a prefix for output files.\n", args->out_fp);
	}

	if (args->doAMOVA != 3 && args->doTest == 1)
	{
		fprintf(stderr, "\n[ERROR]\t-doTest 1 requires -doAMOVA 3; will exit!\n");
		free(args);
		return 0;
	}

	// TODO formatthese text [INFO] [ERROR] etc

	if (args->in_mtd_fn == NULL)
	{
		if (args->doAMOVA != -1)
		{
			fprintf(stderr, "\n[ERROR]\tMust supply -m <metadata_file>; will exit!\n");
			free(args);
			return 0;
		}
	}
	else
	{
		if (args->whichCol == 1)
		{
			fprintf(stderr, "\n[ERROR](-mCol 1)\tColumn index 1 was chosen. First column should contain individual IDs instead; will exit!\n");
			free(args);
			return 0;
		}
		else if (args->whichCol > 1)
		{
			fprintf(stderr, "\n\t-> -mCol is set to %d, will use column %d in metadata file %s as stratification key.\n", args->whichCol, args->whichCol, args->in_mtd_fn);
		}
		else if (args->whichCol == -1)
		{
			args->whichCol = 2;
			fprintf(stderr, "\n\t-> -mCol is not defined, will use column %d in metadata file %s as stratification key.\n", args->whichCol, args->in_mtd_fn);
		}
	}

	// TODO exit(1) or return?
	// maybe dont call these error
	if (args->doDist == -1)
	{
		fprintf(stderr, "\n[ERROR]\tMust supply -doDist <distance_method>; will exit!\n");
		free(args);
		return 0;
	}
	else if (args->doDist == 0)
	{
		fprintf(stderr, "\n\t-> -doDist is set to 0, will use Sij similarity index as distance measure.\n");
	}
	else if (args->doDist == 1)
	{
		fprintf(stderr, "\n\t-> -doDist is set to 1, will use Dij (1-Sij) dissimilarity index as distance measure.\n");
	}
	else if (args->doDist == 2)
	{
		fprintf(stderr, "\n\t-> -doDist is set to 2, will use Fij as distance measure.\n");
	}
	else
	{
		fprintf(stderr, "\n[ERROR]\t-doDist %d is not available; will exit!\n", args->doDist);
		free(args);
		return 0;
	}

	if (args->sqDist == 1)
	{
		fprintf(stderr, "\n\t-> -sqDist is set to 1, will use squared distance measure (dist_ij^2).\n");
	}
	else
	{
		fprintf(stderr, "\n\t-> -sqDist is set to 0, will use absolute value of distance measure (|dist_ij|).\n");
	}

	if (args->doAMOVA == 1)
	{

		if (args->doInd == 1)
		{
			if (args->ind1 == -1)
			{
				fprintf(stderr, "[ERROR]\tMust supply -ind1 while using -doInd 1 \n");
				free(args);
				return 0;
			}
			if (args->ind2 == -1)
			{
				fprintf(stderr, "[ERROR]\tMust supply -ind2 while using -doInd 1 \n");
				free(args);
				return 0;
			}
			if (args->ind1 == args->ind2)
			{
				fprintf(stderr, "[ERROR]\tInd ids must be different while using -doInd 1 \n");
				free(args);
				return 0;
			}
		}

		fprintf(stderr, "\n\t-> -doAMOVA 1; will use 10 genotype likelihoods from GL tag.\n");
	}
	else if (args->doAMOVA == 2)
	{
		fprintf(stderr, "\n\t-> -doAMOVA 2; will use genotypes from GT tag.\n");
	}
	else if (args->doAMOVA == 3)
	{
		fprintf(stderr, "\n\t-> -doAMOVA 2; will do both 1 and 2.\n");
	}
	else if (args->doAMOVA == -1)
	{
		fprintf(stderr, "\n\t-> -doAMOVA -1; will not run AMOVA\n");
	}
	else
	{
		fprintf(stderr, "\n[ERROR]\tMust supply a value for -doAMOVA; will exit!\n");
		free(args);
		return 0;
	}

	return args;
}

paramStruct *paramStruct_init(argStruct *args)
{

	paramStruct *pars = new paramStruct;

	// number of sites non skipped for all individuals
	// nSites may not be !=totSites if minInd is set
	// or if a site is missing for all inds
	pars->nSites = 0;

	// total number of sites processed
	pars->totSites = 0;

	pars->LUT_indPair_idx = NULL;
	pars->n_ind_cmb = 0;

	pars->nInd = 0;

	pars->DATETIME = NULL;

	pars->pos = NULL;

	pars->major = NULL;
	pars->minor = NULL;
	pars->ref = NULL;
	pars->anc = NULL;
	pars->der = NULL;

	return pars;
}

void paramStruct_destroy(paramStruct *pars)
{

	delete[] pars->pos;

	if (pars->major)
	{
		delete[] pars->major;
		pars->major = NULL;
	}
	if (pars->minor)
	{
		delete[] pars->minor;
		pars->minor = NULL;
	}
	delete[] pars->ref;
	delete[] pars->anc;
	delete[] pars->der;

	delete pars->DATETIME;
	pars->DATETIME = NULL;

	for (int i = 0; i < pars->nInd; i++)
	{
		free(pars->LUT_indPair_idx[i]);
		pars->LUT_indPair_idx[i] = NULL;
	}
	free(pars->LUT_indPair_idx);

	delete pars;
}
