/*
 *
 * [Parameters]
 *
 * idea from angsd funkyPars
 * at https://github.com/ANGSD/angsd/blob/master/shared.cpp
 *
 */

#include "shared.h"
#include "math_utils.h"
#include <stdlib.h>

#ifndef FREAD_BUF_SIZE
#define FREAD_BUF_SIZE 4096
#endif

// maximum number of tokens in amova formula string
#ifndef MAX_FORMULA_TOKENS
#define MAX_FORMULA_TOKENS 10
#endif


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


/// @brief get current time
/// @return time as char*
char *get_time()
{
	time_t current_time;
	struct tm *local_time;
	current_time = time(NULL);
	local_time = localtime(&current_time);
	return (asctime(local_time));
}

/// @brief get file handle fp
/// @param fname file name
/// @param mode file open mode
/// @return file *fp
FILE *IO::getFile(const char *fname, const char *mode)
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

/// @brief set file name from prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return filename ie combination of prefix and suffix
char *IO::setFileName(const char *a, const char *b)
{
	char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
	strcpy(c, a);
	strcat(c, b);
	// fprintf(stderr,"\t-> Opening output file for writing: %s\n",c);
	return c;
}

/// @brief open file for writing
/// @param c name of file
/// @return  file *fp
FILE *IO::openFile(char *c)
{
	fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
	FILE *fp = getFile(c, "w");
	return fp;
}


/// @brief open file for writing using given prefix and suffix
/// @param a prefix
/// @param b suffix
/// @return file *fp
FILE *IO::openFile(const char *a, const char *b)
{
	char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
	strcpy(c, a);
	strcat(c, b);
	fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
	FILE *fp = getFile(c, "w");
	free(c);
	return fp;
}


/// @brief getFirstLine get first line of a file
/// @param ff pointer to file
/// @return char* line
char* IO::readFile::getFirstLine(char* in_fn)
{
	FILE* ff = IO::getFile(in_fn, "r");

	size_t buf_size = 1024;
	char* line = (char*) malloc(buf_size);
	// char line[buf_size];

	// char* line = (char*) malloc(buf_size);
	ASSERT(line != NULL);
	
	ASSERT(fgets(line,1024,ff)!=NULL);
	while (fgets(line, 1024, ff) != NULL)
	{
		//check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			//line was fully read
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			//line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char*) realloc(line, buf_size);
			// char *new_line = (char*) realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line=new_line;
		}

	}
	// fprintf(stderr, "\t-> First line of file: %s\n", line);


	// fprintf(stderr, "\t\n-> First line of file: %s\n---\n", line);
	fclose(ff);
	return strdup(line);

}

/// @brief count_nColumns count number of columns in a line
/// @param line pointer to line char
/// @param delims delimiters
/// @return integer number of columns
int IO::inspectFile::count_nColumns(char *line, const char *delims)
{

	int count=0;
	const char *p = line;
	while (*p)
	{
		if (*p == *delims)
			count++;
		p++;
	}
	return count+1;
}

/// @brief count_nRows count number of rows in a file
/// @param ff 
/// @return integer n number of rows
int IO::inspectFile::count_nRows(char* fn, int HAS_HEADER)
{
	FILE* ff = IO::getFile(fn, "r");
	
    char buf[FREAD_BUF_SIZE];
    int n = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, FREAD_BUF_SIZE, ff);
        if (ferror(ff))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n') n++;

        if (feof(ff))
            break;
    }

	if (HAS_HEADER==1) n--;
	fclose(ff);
    return n;
}

int IO::validateFile::Metadata(char* in_mtd_fn, int nInds, int* keyCols, DATA::formulaStruct* fos, const char* delims, int HAS_HEADER)
{
	int nRows = IO::inspectFile::count_nRows(in_mtd_fn, HAS_HEADER);
	if (nRows == -1)
	{
		fprintf(stderr, "\n[ERROR]\tCould not count number of rows in Metadata file; will exit!\n\n");
		exit(1);
	}
	if (nRows == 0)
	{
		fprintf(stderr, "\n[ERROR]\tMetadata file is empty; will exit!\n\n");
		exit(1);
	}
	if (nRows == 1)
	{
		fprintf(stderr, "\n[ERROR]\tMetadata file contains only one row; will exit!\n\n");
		exit(1);
	}
	if (nRows != nInds)
	{
		fprintf(stderr, "\n[ERROR]\tNumber of rows in Metadata file (%d) does not match number of individuals (%d); will exit!\n\n", nRows, nInds);
		exit(1);
	}




	//compare number of tokens in formulaStruct to number of columns in Metadata file

	char* firstLine = IO::readFile::getFirstLine(in_mtd_fn);
	int nCols = IO::inspectFile::count_nColumns(firstLine, delims);
	fprintf(stderr, "\n\t-> Number of columns in input Metadata file: %d\n", nCols);

	if (nCols < fos->nTokens)
	{
		fprintf(stderr, "\n[ERROR]\tNumber of columns in Metadata file (%d) is less than number of tokens in formula (%d); will exit!\n\n", nCols, fos->nTokens);
		exit(1);
	}

	// keyCols

	// int nCols = IO::inspectFile::count_nColumns(IO::inspectFile::getLine(ff, 1), "\t ");

	return 0;


}

/// @brief read SFS file
/// @param in_sfs_ff input sfs file ff
/// @param delims delimiters
/// @param SAMPLES samplesStruct samples
/// @return ???
int IO::readFile::SFS(FILE *in_sfs_ff, const char *delims, DATA::samplesStruct *SAMPLES)
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


/// @brief read Metadata file
/// @param MTD metadataStruct Metadata
/// @param in_mtd_fn input Metadata file filename
/// @param in_mtd_ff input Metadata file ff
/// @param keyCols pointer to set of which columns to use as keys
/// @param delims delimiters
/// @param SAMPLES samplesStruct samples
/// @param FORMULA formulaStruct formula
/// @param HAS_HEADER metadata file has header, yes=1, no=0
/// @return ???
int IO::readFile::Metadata(DATA::metadataStruct *MTD, char* in_mtd_fn, FILE *in_mtd_ff, int* keyCols, const char *delims, DATA::samplesStruct *SAMPLES, DATA::formulaStruct *FORMULA, int HAS_HEADER)
{

	fprintf(stderr, "\n");
	fprintf(stderr, "--------------------------------------------------");
	fprintf(stderr, "\n");

	// int nCols = 0;
	int checkCol = 1;

	// int* keyCols = (int*)malloc(sizeof(int) * FORMULA->nTokens);
	fprintf(stderr,"\n\nnTokens: %d\n", FORMULA->nTokens);
	fprintf(stderr, "\n formula tokens: ");
	for (int i = 0; i < FORMULA->nTokens; i++)
	{
		fprintf(stderr, "%s ", FORMULA->formulaTokens[i]);
	}
	fprintf(stderr, "\n");

	ASSERT(validateFile::Metadata(in_mtd_fn, SAMPLES->nSamples, keyCols, FORMULA, delims, HAS_HEADER) == 0);
	// "Metadata file is not valid");

	// TODO map is probably better due to sorting issues. we cannot have all strata sorted
	char mt_buf[1024];

	//skip the first line
	if (HAS_HEADER)
	{
		fgets(mt_buf, 1024, in_mtd_ff);
	}


	while (fgets(mt_buf, 1024, in_mtd_ff))
	{



		char *tok = strtok(mt_buf, delims);
		char *group_id = tok;

//TODO deprecated
		int whichCol=-1;
		for (int coli = 0; coli < whichCol - 1; coli++)
		{
			tok = strtok(NULL, "\t \n");
			group_id = tok;
		}

		// increase the size of Strata
		if (MTD->nStrata > (int)MTD->_strataArr)
		{
			fprintf(stderr, "->->->increase the size of Strata S[]!! Found %d MTD->nStrata and %d MTD->_strataArr\n", (int)MTD->nStrata, (int)MTD->_strataArr);
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

// IO::print::Array
/// @brief Print elements of an array to a file (or stream)
/// @param fp file to print to
/// @param arr the array to print
/// @param N the number of rows, dimension 1
/// @param M the number of columns, dimension 2
/// @param sep the character to use as a separator
/// @example IO::print::Array(stdout, myArray, DOUBLE, 3, 3, ',');
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

/// @brief file_exists - check if file exists
/// @param in_fn input filename
/// @return 1 if file exists; 0 otherwise
/// @credit angsd/aio.cpp
int file_exists(const char *in_fn)
{
	struct stat buffer;
	return (stat(in_fn, &buffer) == 0);
}

/// @brief usage - print usage
/// @param fp pointer to the file to print to
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

	// fprintf(fp, "\n");
	// fprintf(fp, "Usage: ngsAMOVA [options] -in <vcf file> -out <output file>\n");
	// fprintf(fp, "\n");
	// fprintf(fp, "Options:\n");
	// fprintf(fp, "  -in <vcf file>		: input vcf file\n");
	// fprintf(fp, "  -out <output file>		: output file\n");
	// fprintf(fp, "  -doAMOVA <0/1>		: do AMOVA (default: 1)\n");
	// fprintf(fp, "  -doTest <0/1>		: do EM test (default: 1)\n");
	// fprintf(fp, "  -printMatrix <0/1>		: print distance matrix (default: 0)\n");
	// fprintf(fp, "  -printSFS <0/1>		: print SFS (default: 0)\n");


	//
	//
	// exit(0);
}



/// @brief argStruct_init - initialize the argStruct structure
/// @return pointer to the argStruct structure
argStruct *argStruct_init()
{

	argStruct *args = (argStruct *)calloc(1, sizeof(argStruct));

	args->in_fn = NULL;
	args->in_sfs_fn = NULL;
	args->in_mtd_fn = NULL;
	args->out_fp = NULL;

	args->formula = NULL;
	args->keyCols = NULL;

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
	args->isTest = 0;
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
		else if (strcasecmp("-af", arv) == 0){
			args->formula = strdup(val);
		}
		else if (strcasecmp("--formula", arv) == 0){
			args->formula = strdup(val);
		}
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
		else if (strcasecmp("-isTest", arv) == 0)
			args->isTest = atoi(val);
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

	if (args->isTest == 1)
	{
		fprintf(stderr,"Test mode ON\n");
	}

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
			fprintf(stderr, "\n[ERROR]\tMust supply -m <Metadata_file>; will exit!\n");
			free(args);
			return 0;
		}
	}
	else
	{
		// if (args->keyCols == 1)
		// {
		// 	fprintf(stderr, "\n[ERROR](-mCol 1)\tColumn index 1 was chosen. First column should contain individual IDs instead; will exit!\n");
		// 	free(args);
		// 	return 0;
		// }
		// else if (args->keyCols > 1)
		// {
		// 	fprintf(stderr, "\n\t-> -mCol is set to %d, will use column %d in Metadata file %s as stratification key.\n", args->keyCols, args->keyCols, args->in_mtd_fn);
		// }
		// else if (args->keyCols == -1)
		// {
		// 	args->keyCols = 2;
		// 	fprintf(stderr, "\n\t-> -mCol is not defined, will use column %d in Metadata file %s as stratification key.\n", args->keyCols, args->in_mtd_fn);
		// }
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
		exit(1);
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
		exit(1);
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


/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
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

/// @brief paramStruct_destroy free memory of paramStruct
/// @param pars pointer to paramStruct
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

DATA::formulaStruct *DATA::formulaStruct_init(const char *formula)
{

	formulaStruct *fos = (formulaStruct *)calloc(1, sizeof(formulaStruct));
	
	
	// char** formulaTokens[MAX_FORMULA_TOKENS];
	// formula = 'Samples ~ Continents/Regions/Populations'

	fos->formulaTokens = (char**)malloc(MAX_FORMULA_TOKENS * sizeof(char*));
	fos->nTokens = 0;

	fos->formula = strdup(formula);
	// pointer to the first character of formula string
	const char *p = formula;

	int nTilde=0;

	// skip until ~ tilde or end of string
	while (*p != '\0' && *p != '~')
	{
		if (*p == '/') {
			fprintf(stderr, "\n[ERROR]\tformula %s is not valid: Found '/' before '~'. Will exit!\n", formula);
			exit(1);
		}
		if (*p == '~') {
			nTilde++;
		}
		p++;
	}

	if(*p == '\0'){
		fprintf(stderr, "\n[ERROR]\tformula %s is not valid; will exit!\n", formula);
		exit(1);
	}

	fos->formulaTokens[fos->nTokens] = strndup(formula, p - formula);
	fprintf(stderr, "\n\n\n-------\nntokens: %d\n", fos->nTokens);
	fos->nTokens++;
	p++;
	
	while (*p != '\0')
	{
		//start of current token
		const char* start = p;
		
		//skip until '/' is encountered or end of string
		while(*p != '\0' && *p != '/')
		{
			p++;
		}
		//end of string
		if(*p == '\0'){
			fos->formulaTokens[fos->nTokens]= strndup(start, p - start);
			fos->nTokens++;
			break;
		}
		if (*p == '~') {
			nTilde++;
			ASSERT(nTilde!=1); //should never happen
			if(nTilde > 1){
				fprintf(stderr, "\n[ERROR]\tformula %s is not valid: Found more than one '~'. Will exit!\n", formula);
			}
		}
		
		//or it was a '/'
		fos->formulaTokens[fos->nTokens]= strndup(start, p - start);
		fos->nTokens++;
		p++;
	}
	// set the last element to NULL to indicate end of array
	fos->formulaTokens[fos->nTokens] = NULL;




	// int *keyCols = (int *)malloc(fos->nTokens * sizeof(int));

	


	return fos;

}

// void IO::readArgs::printArgs(argStruct *args)
// {

// 	fprintf(stderr, "\n\t-> -in %s", args->inputfile);
// }


/// @brief print_M_PWD print matrix of pairwise distances
/// @param TYPE type of analysis
/// @param out_dm_fs output file
/// @param n_ind_cmb number of individual combinations
/// @param M_PWD matrix of pairwise distances
void IO::print::M_PWD(const char *TYPE, IO::outputStruct *out_dm_fs, int n_ind_cmb, double *M_PWD)
{

	fprintf(out_dm_fs->ff, "%s,", TYPE);
	for (int px = 0; px < n_ind_cmb; px++)
	{
		fprintf(out_dm_fs->ff, "%f", M_PWD[px]);
		if (px != n_ind_cmb - 1)
		{
			fprintf(out_dm_fs->ff, ",");
		}
		else
		{
			fprintf(out_dm_fs->ff, "\n");
		}
	}
}



/// @param sample1 name of sample 1
/// @param sample2 name of sample 2
void IO::print::Sfs(const char* TYPE, IO::outputStruct *out_sfs_fs, DATA::pairStruct *pair, argStruct *args,const char* sample1, const char* sample2)
{


		fprintf(out_sfs_fs->ff, "%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f",
				TYPE,
				sample1,
				sample2,
				pair->snSites * pair->SFS[0], pair->snSites * pair->SFS[1], pair->snSites * pair->SFS[2],
				pair->snSites * pair->SFS[3], pair->snSites * pair->SFS[4], pair->snSites * pair->SFS[5],
				pair->snSites * pair->SFS[6], pair->snSites * pair->SFS[7], pair->snSites * pair->SFS[8]);

		fprintf(out_sfs_fs->ff, ",%d,%ld,%e,%e", pair->n_em_iter, pair->snSites, pair->d, args->tole);

		if(args->doDist == 1 && args->sqDist == 1){
			fprintf(out_sfs_fs->ff, ",%f,%f", (double) (1.0 - (double) MATH::EST::Sij(pair->SFS) ), SQUARE((double) (1.0 - (double) MATH::EST::Sij(pair->SFS) )));
		}else{
			exit(1);
		}
		fprintf(out_sfs_fs->ff,"\n");

}


/// @brief print_SFS_GT print SFS_GT3
/// @param TYPE type of analysis
/// @param out_sfs_fs output file
/// @param args pointer to argStruct
/// @param SFS_GT3 matrix of 3 GT SFS for pair (int **SFS_GT3[pidx])
/// @param snSites (shared) number of sites
void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char* sample1, const char* sample2)
{

	fprintf(out_sfs_fs->ff, "%s,%s,%s,", TYPE, sample1, sample2);
	fprintf(out_sfs_fs->ff, "%d,%d,%d,%d,%d,%d,%d,%d,%d",
			SFS_GT3[0], SFS_GT3[1], SFS_GT3[2],
			SFS_GT3[3], SFS_GT3[4], SFS_GT3[5],
			SFS_GT3[6], SFS_GT3[7], SFS_GT3[8]);
	
	fprintf(out_sfs_fs->ff, ",%s,%ld,%s,%s", TYPE, snSites, TYPE, TYPE);

	if(args->doDist == 1 && args->sqDist == 1){
			fprintf(out_sfs_fs->ff, ",%f", SQUARE((double) (1.0 - (double) MATH::EST::Sij(SFS_GT3, snSites) )));
	}else{
		exit(1);
	}
	fprintf(out_sfs_fs->ff, "\n");
}
