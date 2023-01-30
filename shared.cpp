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

// TODO check this
using size_t = decltype(sizeof(int));


int find_n_given_nC2(int nC2_res){
	int n=0;
	while(nChoose2[n]<nC2_res){
		n++;
	}
	if(nChoose2[n]!=nC2_res){
		fprintf(stderr,"[%s:%s()]\t->Error: nC2_res:%d not found in nChoose2[]\n",__FILE__,__FUNCTION__,nC2_res);
		exit(1);
	}
	return n;
}



//TODO
// extract digits using bit masking
// int extractDigits(int num, int digits)
// {
// 	uint32_t x = 100101;
// 	uint32_t y = 100102;
// 	for(int i=0;i<sizeof(uint32_t);i++){
// 		fprintf(stderr,"\n\n\n\n -> x: %d y: %d\n",(x>>i)&0xF,(y>>i)&0xF);
// 	}
// 	exit(0);
// }



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
/// @param fn file name
/// @return char* line
char *IO::readFile::getFirstLine(char *fn)
{
	FILE *fp = IO::getFile(fn, "r");

	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);

	ASSERT(line != NULL);

	// ASSERT(fgets(line, 1024, fp) != NULL);
	while (fgets(line, 1024, fp) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			// line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char *)realloc(line, buf_size);
			// char *new_line = (char*) realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line = new_line;
		}
	}
	// fprintf(stderr, "\t-> First line of file: %s\n", line);
	fclose(fp);
	// TODO check this return
	return strdup(line);
}

/// @brief getFirstLine get first line of a file
/// @param fp pointer to file
/// @return char* line
char *IO::readFile::getFirstLine(FILE *fp)
{

	ASSERT(fseek(fp, 0, SEEK_SET) == 0);
	size_t buf_size = FGETS_BUF_SIZE;
	char *line = (char *)malloc(buf_size);

	ASSERT(line != NULL);

	while (fgets(line, 1024, fp) != NULL)
	{
		// check if the line was fully read
		size_t line_len = strlen(line);
		if (line[line_len - 1] == '\n')
		{
			// line was fully read
			break;
		}
		else
		{
			fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
			// line was not fully read
			buf_size *= 2;

			char *new_line = new char[buf_size];
			new_line = (char *)realloc(line, buf_size);
			// char *new_line = (char*) realloc(line, buf_size);
			ASSERT(new_line != NULL);
			line = new_line;
		}
	}
	// fprintf(stderr, "\t-> First line of file: %s\n", line);
	return strdup(line);
}

/// @brief count_nColumns count number of columns in a line
/// @param line pointer to line char
/// @param delims delimiters
/// @return integer number of columns
int IO::inspectFile::count_nColumns(char *line, const char *delims)
{

	int count = 0;
	const char *p = line;
	while (*p)
	{
		if (*p == *delims)
			count++;
		p++;
	}
	return count + 1;
}

/// @brief count_nRows count number of rows in a file
/// @param fn file name
/// @param HAS_COLNAMES 1 if file has header
/// @return integer n number of rows
int IO::inspectFile::count_nRows(char *fn, int HAS_COLNAMES)
{
	FILE *fp = IO::getFile(fn, "r");

	char buf[FREAD_BUF_SIZE];
	int n = 0;
	for (;;)
	{
		size_t res = fread(buf, 1, FREAD_BUF_SIZE, fp);
		ASSERT(ferror(fp) == 0);

		size_t i;
		for (i = 0; i < res; i++)
			if (buf[i] == '\n')
				n++;

		if (feof(fp))
			break;
	}

	if (HAS_COLNAMES == 1)
		n--;

	fclose(fp);

	return n;
}

/// @brief count_nRows count number of rows in a file
/// @param fp pointer to file
/// @param HAS_COLNAMES 1 if file has header
/// @return integer n number of rows
int IO::inspectFile::count_nRows(FILE *fp, int HAS_COLNAMES)
{
	// return to the beginning of the file
	ASSERT(fseek(fp, 0, SEEK_SET) == 0);

	char buf[FREAD_BUF_SIZE];
	int n = 0;
	for (;;)
	{
		size_t res = fread(buf, 1, FREAD_BUF_SIZE, fp);
		ASSERT(ferror(fp) == 0);

		size_t i;
		for (i = 0; i < res; i++)
			if (buf[i] == '\n')
				n++;

		if (feof(fp))
			break;
	}

	if (HAS_COLNAMES == 1)
		n--;
	return n;
}

/// @brief IO::validateFile::Metadata validate metadata file (input=BCF)
/// @param in_mtd_fp pointer to metadata file
/// @param nInds number of individuals
/// @param keyCols key columns
/// @param FORMULA formula
/// @param delims delimiters
/// @param HAS_COLNAMES 1 if file has header
/// @return number of columns if successful, exits with error otherwise
int IO::validateFile::Metadata(FILE *in_mtd_fp, int nInds, int *keyCols,
							   DATA::formulaStruct *FORMULA, const char *delims, int HAS_COLNAMES)
{

	ASSERT(fseek(in_mtd_fp, 0, SEEK_SET) == 0);
	int nRows = IO::inspectFile::count_nRows(in_mtd_fp, HAS_COLNAMES);
	if (nRows == -1)
	{
		fprintf(stderr, "\n[ERROR]\tCould not count number of rows in Metadata file.\n\n");
		exit(1);
	}
	if (nRows == 0)
	{
		fprintf(stderr, "\n[ERROR]\tMetadata file is empty.\n\n");
		exit(1);
	}
	if (nRows == 1)
	{
		fprintf(stderr, "\n[ERROR]\tMetadata file contains only one row.\n\n");
		exit(1);
	}
	if (nRows != nInds)
	{
		fprintf(stderr, "\n[ERROR]\tNumber of rows in Metadata file (%d) does not match number of individuals (%d).\n\n", nRows, nInds);
		exit(1);
	}

	// compare number of tokens in formulaStruct to number of columns in Metadata file

	char *firstLine = IO::readFile::getFirstLine(in_mtd_fp);
	int nCols = IO::inspectFile::count_nColumns(firstLine, delims);
	fprintf(stderr, "\n\t-> Number of columns in input Metadata file: %d\n", nCols);

	if (FORMULA != NULL)
	{
		if (nCols < FORMULA->nTokens)
		{
			fprintf(stderr, "\n[ERROR]\tNumber of columns in Metadata file (%d) is less than number of tokens in formula (%d).\n\n", nCols, FORMULA->nTokens);
			exit(1);
		}
	}

	// int nCols = IO::inspectFile::count_nColumns(IO::inspectFile::getLine(ff, 1), "\t ");
	return nCols;
}

/// @brief read SFS file
/// @param in_sfs_fp input sfs file ff
/// @param delims delimiters
/// @param sampleSt sampleStruct samples
/// @return ???
int IO::readFile::SFS(FILE *in_sfs_fp, const char *delims, DATA::sampleStruct *sampleSt)
{

	char sfs_buf[FGETS_BUF_SIZE];
	while (fgets(sfs_buf, FGETS_BUF_SIZE, in_sfs_fp))
	{

		char *tok = strtok(sfs_buf, delims);
		char *col = tok;

		for (int coli = 0; coli < 9 - 1; coli++)
		{
			tok = strtok(NULL, delims);
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

/// @brief read distance matrix file
/// @param in_dm_fp input distance matrix file
/// @param pars paramStruct parameters
/// @return distance matrix double*
DATA::distanceMatrixStruct *DATA::distanceMatrixStruct_read(FILE *in_dm_fp, paramStruct *pars, argStruct *args)
{

	char dm_buf[FGETS_BUF_SIZE];

	int dm_vals_size = 1225;
	double *dm_vals = (double *)malloc((dm_vals_size) * sizeof(double));

	// first value is type string indicating analysis type
	// TODO exclude this from distance matrix
	// int n_vals = -1;
	int n_vals = 0;

	while (fgets(dm_buf, FGETS_BUF_SIZE, in_dm_fp))
	{
		char *tok = strtok(dm_buf, ",");
		while (tok != NULL)
		{
			if (n_vals > dm_vals_size)
			{
				dm_vals_size = dm_vals_size * 2;
				dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
			}

			dm_vals[n_vals] = atof(tok);
			n_vals++;
			tok = strtok(NULL, ", \n");
		}
	}
	pars->nIndCmb = n_vals;
	fprintf(stderr, "\t-> Number of values in distance matrix: %d. (i.e. number of unique individual pairs)\n", n_vals);
	pars->nInd = find_n_given_nC2(n_vals);
	fprintf(stderr, "\n\t-> Number of individuals in distance matrix: %d\n", pars->nInd);
	pars->init_LUTs();
	set_LUT_inds2idx_2way(pars->nInd, pars->nIndCmb, pars->LUT_inds2idx, pars->LUT_idx2inds);

	distanceMatrixStruct *dMS = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->do_square_distance);
	dMS->isSquared = args->do_square_distance;

	if (args->do_square_distance == 1)
	{
		for (int i = 0; i < n_vals; i++)
		{
			dMS->M[i] = SQUARE(dm_vals[i]);
		}
	}
	else
	{
		for (int i = 0; i < n_vals; i++)
		{
			dMS->M[i] = dm_vals[i];
		}
	}

	check_consistency_args_pars(args, pars);

	FREE(dm_vals);


	return dMS;
}


// if formula is not defined, hierarchical levels must defined in this order:
// individual, highest_level, ..., lowest_level
// e.g.
// individual, region, population, subpopulation
// which translates to the formula:
// individual ~ region / population / subpopulation
DATA::metadataStruct *DATA::metadataStruct_get(FILE *in_mtd_fp, DATA::sampleStruct *sampleSt,
											   DATA::formulaStruct *FORMULA, int has_colnames,
											   paramStruct *pars)
{

	char *mtd_buf = new char[FGETS_BUF_SIZE];

	size_t nLevels = 0;

	// go to beginning of file
	ASSERT(fseek(in_mtd_fp, 0, SEEK_SET) == 0);

	if (FORMULA == NULL)
	{

		if (has_colnames != 1){
			fprintf(stderr, "\n[ERROR]\t If no formula is provided, the metadata file must have a header line with the column names.\n");
			exit(1);
		}
			
		// if no formula is provided, then we assume that columns are ordered based on hierarchical structure
		fprintf(stderr, "\n[INFO]\t-> Formula is not set, will assume that columns are ordered based on hierarchical structure and will use all columns.\n");

		// names of levels in the hierarchy (e.g. region, population, subpopulation)
		char **levelNames = new char *[MAX_N_AMOVA_LEVELS];

		// skip first line
		ASSERT(fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp) != NULL);

		// first column contains the keyword corresponding to the individual id
		char *hdrtok = strtok(mtd_buf, METADATA_DELIMS);
		ASSERT(hdrtok != NULL);

		// exclude Individual column
		nLevels = -1;

		do
		{
			fprintf(stderr, "\n[INFO]\t-> Found hierarchical level: %s\n", hdrtok);
			nLevels++;
			levelNames[nLevels] = new char[strlen(hdrtok) + 1];
			ASSERT(strncpy(levelNames[nLevels], hdrtok, strlen(hdrtok) + 1) != NULL);
			hdrtok = strtok(NULL, METADATA_DELIMS);
		} while (hdrtok != NULL);

		if (nLevels > MAX_N_AMOVA_LEVELS)
		{
			fprintf(stderr, "\n[ERROR]\t-> Number of levels in metadata file (%ld) exceeds the maximum number of levels allowed (%d).\n", nLevels, MAX_N_AMOVA_LEVELS);
			exit(1);
		}

		fprintf(stderr, "\n[INFO]\t-> Number of levels in metadata file: %ld\n", nLevels);


		metadataStruct *mS = new metadataStruct(nLevels, levelNames);

		// collect the strata_i indexes for each hierarchical level until the lowest level
		size_t key = 0;

		// sample index 
		int sidx= 0;

		// loop through the remaining lines; one line per individual
		while (fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp))
		{

			mS->nIndMetadata++;

			// strata key
			// e.g. 4 levels: continent, region, population, subpopulation
			// key init 1e ((4-1)*2) -> 1e6 -> 1000000
			// key 1000310 represents one subpopulation
			//     1------
			//     -00 (belongs to the continent with index 0, i.e. continent1)
			//	   ---03 (belongs to the region with index 3, i.e. region4)
			//     -----10 (belongs to the population with index 10, i.e. population11)
			// key = (size_t) 1 * (size_t) pow(10, MAXDIG_PER_HLEVEL * nLevels);
			key = mS->initKey();

			// first column = individual id
			char *tok = strtok(mtd_buf, METADATA_DELIMS);



			// ----------------------- INPUT: VCF ----------------------- //
			if (pars->in_ft == IN_VCF)
			{

				// if input file type is vcf, find the index of the individual in the bcf header
				for (sidx = 0; sidx < sampleSt->nSamples + 1; sidx++)
				{

					if (sidx == sampleSt->nSamples)
					{
						fprintf(stderr, "\n\n======\n[ERROR] Sample %s not found in the metadata file. \n\n", tok);
						exit(1);
					}
					else if (strcmp(tok, sampleSt->sampleNames[sidx]) == 0)
					{
						if (pars->verbose == 2) {
							fprintf(stderr, "\n\n======\n[INFO] Found sample (vcf_index=%d,id=%s) in the metadata file. \n\n", sidx, tok);
						}
						// break when found, so that sidx is now the index of the sample in the vcf file
						break;
					}
				}

				// loop through the hierarchical levels in the line
				for (int lvl_i = 0; lvl_i < mS->nLevels; lvl_i++)
				{

					//index of the strata in the current hierarchical level
					int lvl_strata_i=0;

					// first column after individual is the highest hierarchical level
					tok = strtok(NULL, METADATA_DELIMS);


					// if we are reading the first line
					if (mS->hierArr[lvl_i] == NULL)
					{
						// we are in the first line, initialize the hierStruct for this level

						//TODO
						// mS->addHierStruct(lvl_i, tok);
						mS->hierArr[lvl_i] = new hierStruct(tok);
						mS->hierArr[lvl_i]->nIndPerStrata[0]++;


						// EOL end of the first line
						if (lvl_i == mS->nLevels - 1)
						{
							// EOL
							// lowest level (excluding individual) in the first line
							// no check because this is the first line of data we read
							ASSERT(lvl_strata_i == 0);
							key = mS->setKeyDigitAtLevel(key, lvl_i,lvl_strata_i);
							mS->ind2stratakey[sidx] = key;
							continue;
						}

					}

					// EOL we are at the lowest level column
					// associate the individual with the stratakey
					else if (lvl_i == mS->nLevels - 1)
					{

						int strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);

						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);

						mS->ind2stratakey[sidx] = key;


						// at the lowest level; end of the individual's line

					}
					else
					{
						// at a regular column, nothing special about it
						// (not the first line; not the first column; not the last column)
						//
						// e.g. ind1,x1,y1,z1
						// 		ind2,[x1],y2,z3  
						//   		we are at 'x1' of ind2, already allocated in previous individual

						int strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
					}
				}
			}
			// ----------------------- INPUT: DM ----------------------- //
			else if (pars->in_ft == IN_DM)
			{

				sampleSt->addSampleName(sidx, tok);
				// loop through the hierarchical levels in the line
				for (int lvl_i = 0; lvl_i < mS->nLevels; lvl_i++)
				{

					//index of the strata in the current hierarchical level
					int lvl_strata_i=0;

					// first column after individual is the highest hierarchical level
					tok = strtok(NULL, METADATA_DELIMS);


					// if we are reading the first line
					if (mS->hierArr[lvl_i] == NULL)
					{
						// we are in the first line, initialize the hierStruct for this level

						//TODO
						// mS->addHierStruct(lvl_i, tok);
						mS->hierArr[lvl_i] = new hierStruct(tok);
						mS->hierArr[lvl_i]->nIndPerStrata[0]++;


						// EOL end of the first line
						if (lvl_i == mS->nLevels - 1)
						{
							// EOL
							// lowest level (excluding individual) in the first line
							// no check because this is the first line of data we read
							ASSERT(lvl_strata_i == 0);
							key = mS->setKeyDigitAtLevel(key, lvl_i,lvl_strata_i);
							mS->ind2stratakey[sidx] = key;
							sidx++;
							continue;
						}

					}

					// EOL 
					// at the lowest level column; end of the individual's line
					else if (lvl_i == mS->nLevels - 1)
					{

						// associate the individual with the stratakey

						int strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);

						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);

						mS->ind2stratakey[sidx] = key;

						sidx++;

					}
					else
					{
						// at a regular column, nothing special about it
						// (not the first line; not the first column; not the last column)
						//
						// e.g. ind1,x1,y1,z1
						// 		ind2,[x1],y2,z3  
						//   		we are at 'x1' of ind2, already allocated in previous individual

						int strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
					}
				}
			}
			else
			{
				fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
				ASSERT(0 == 1);
			}
		}

		for (size_t i = 0; i < (size_t)mS->nLevels + 1; i++)
		{
			delete[] levelNames[i];
		}
		delete[] levelNames;
		delete[] mtd_buf;


		return (mS);

	}else{
		// formula is defined, i.e. not NULL

		fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
		ASSERT(0 == 1);
	}

	ASSERT(0 == 1);
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
// @param fn input filename
// @return		1 if file exists; 0 otherwise

/// @brief file_exists - check if file exists
/// @param fn input filename
/// @return 1 if file exists; 0 otherwise
/// @credit angsd/aio.cpp
int file_exists(const char *fn)
{
	struct stat buffer;
	return (stat(fn, &buffer) == 0);
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
	// fprintf(fp, "  -doTest <0/1>		: do EM test (default: 
	// fprintf(fp, "  -printMatrix <0/1>		: print distance matrix (default: 0)\n");
	// fprintf(fp, "  -printSFS <0/1>		: print SFS (default: 0)\n");

	//
	//
	exit(0);
}

/// @brief argStruct_init - initialize the argStruct structure
/// @return pointer to the argStruct structure
argStruct *argStruct_init()
{

	argStruct *args = (argStruct *)calloc(1, sizeof(argStruct));

	args->verbose = 1;

	args->in_vcf_fn = NULL;
	args->in_sfs_fn = NULL;
	args->in_dm_fn = NULL;
	args->in_mtd_fn = NULL;
	args->out_fn = NULL;

	args->formula = NULL;
	args->keyCols = NULL;
	args->hasColNames = 1;

	args->command = NULL;
	args->blockSize = 0;

	args->doAMOVA = 0;
	args->doEM = 0;

	args->mThreads = 0;
	args->mEmIter = 1e2;

	args->tole = 1e-5;
	args->doTest = 0;

	args->doDist = -1;
	args->do_square_distance = 1;
	args->printMatrix = 0;

	args->isSim = 0;
	args->isTest = 0;
	args->minInd = -1;

	args->printDev= -1;


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

		if (strcasecmp("-in", arv) == 0)
			args->in_vcf_fn = strdup(val);
		else if (strcasecmp("-i", arv) == 0)
			args->in_vcf_fn = strdup(val);
		else if (strcasecmp("-in_sfs", arv) == 0)
			args->in_sfs_fn = strdup(val);
		else if (strcasecmp("-in_dm", arv) == 0)
			args->in_dm_fn = strdup(val);
		else if (strcasecmp("-m", arv) == 0)
			args->in_mtd_fn = strdup(val);
		else if (strcasecmp("-out", arv) == 0)
			args->out_fn = strdup(val);
		else if (strcasecmp("-o", arv) == 0)
			args->out_fn = strdup(val);

		else if (strcasecmp("-dev", arv) == 0)
			args->printDev = atoi(val);



		// read block size as float and convert to int
		// this is to allow for the use of scientific notation (e.g. 1e6)
		else if (strcasecmp("-bs", arv) == 0)
			args->blockSize = (int) atof(val);
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

	return args;
}

void argStruct_destroy(argStruct *args)
{
	FREE(args->in_vcf_fn);
	FREE(args->in_sfs_fn);
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


/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct *paramStruct_init(argStruct *args)
{

	paramStruct *pars = new paramStruct;

	if (args->in_vcf_fn != NULL)
	{
		pars->in_ft = IN_VCF;
	}
	else if (args->in_dm_fn != NULL)
	{
		pars->in_ft = IN_DM;
	}

	// number of sites non skipped for all individuals
	// nSites may not be !=totSites if minInd is set
	// or if a site is missing for all inds
	pars->nSites = 0;

	// total number of sites processed
	pars->totSites = 0;

	pars->LUT_inds2idx = NULL;
	pars->LUT_idx2inds = NULL;

	pars->nIndCmb = 0;
	pars->nInd = 0;

	pars->DATETIME = NULL;

	pars->nAmovaRuns = 0;


	pars->verbose=0;

	// pars->pos = NULL;
	// pars->major = NULL;
	// pars->minor = NULL;
	// pars->ref = NULL;
	// pars->anc = NULL;
	// pars->der = NULL;

	return pars;
}

/// @brief paramStruct_destroy free memory of paramStruct
/// @param pars pointer to paramStruct
void paramStruct_destroy(paramStruct *pars)
{

	// delete[] pars->pos;

	// if (pars->major)
	// {
	// 	delete[] pars->major;
	// 	pars->major = NULL;
	// }
	// if (pars->minor)
	// {
	// 	delete[] pars->minor;
	// 	pars->minor = NULL;
	// }
	// delete[] pars->ref;
	// delete[] pars->anc;
	// delete[] pars->der;

	FREE(pars->DATETIME);

	for (int i = 0; i < pars->nInd; i++)
	{
		FREE(pars->LUT_inds2idx[i]);
	}

    for (int i = 0; i < pars->nIndCmb; i++)
    {
        FREE(pars->LUT_idx2inds[i]);
    }

	FREE(pars->LUT_inds2idx);
	FREE(pars->LUT_idx2inds);

	delete pars;
}

/// @brief formulaStruct_get initialize the formulaStruct
/// @param formula formula string
/// @return pointer to formulaStruct
/// @example formula = 'Samples ~ Continents/Regions/Populations'
DATA::formulaStruct *DATA::formulaStruct_get(const char *formula)
{

	formulaStruct *fos = new formulaStruct;

	// formulaStruct *fos = (formulaStruct *)calloc(1, sizeof(formulaStruct));

	fos->formulaTokens = (char **)malloc(MAX_FORMULA_TOKENS * sizeof(char *));
	fos->nTokens = 0;

	fos->formula = strdup(formula);

	// pointer to the first character of formula string
	const char *p = formula;

	int nTilde = 0;

	// skip until ~ tilde or end of string
	while (*p != '\0' && *p != '~')
	{
		if (*p == '/')
		{
			fprintf(stderr, "\n[ERROR]\tformula %s is not valid: Found '/' before '~'. \n", formula);
			exit(1);
		}
		if (*p == '~')
		{
			nTilde++;
		}
		p++;
	}

	if (*p == '\0')
	{
		fprintf(stderr, "\n[ERROR]\tformula %s is not valid.\n", formula);
		exit(1);
	}

	fos->formulaTokens[fos->nTokens] = strndup(formula, p - formula);
	// fprintf(stderr, "\n\n\n-------\nntokens: %d\n", fos->nTokens);
	fos->nTokens++;
	p++;

	while (*p != '\0')
	{
		// start of current token
		const char *start = p;

		// skip until '/' is encountered or end of string
		while (*p != '\0' && *p != '/')
		{
			p++;
		}
		// end of string
		if (*p == '\0')
		{
			fos->formulaTokens[fos->nTokens] = strndup(start, p - start);
			fos->nTokens++;
			break;
		}
		if (*p == '~')
		{
			nTilde++;
			ASSERT(nTilde != 1); // should never happen
			if (nTilde > 1)
			{
				fprintf(stderr, "\n[ERROR]\tformula %s is not valid: Found more than one '~'. \n", formula);
			}
		}
		// or it was a '/'
		fos->formulaTokens[fos->nTokens] = strndup(start, p - start);

		fos->nTokens++;
		p++;
	}

	// int *keyCols = (int *)malloc(fos->nTokens * sizeof(int));

	return fos;
}

// void IO::readArgs::printArgs(argStruct *args)
// {

// 	fprintf(stderr, "\n\t-> -in %s", args->inputfile);
// }

/// @param sample1 name of sample 1
/// @param sample2 name of sample 2
void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, DATA::pairStruct *pair, argStruct *args, const char *sample1, const char *sample2)
{

	fprintf(out_sfs_fs->fp, "%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f",
			TYPE,
			sample1,
			sample2,
			pair->snSites * pair->SFS[0], pair->snSites * pair->SFS[1], pair->snSites * pair->SFS[2],
			pair->snSites * pair->SFS[3], pair->snSites * pair->SFS[4], pair->snSites * pair->SFS[5],
			pair->snSites * pair->SFS[6], pair->snSites * pair->SFS[7], pair->snSites * pair->SFS[8]);

	fprintf(out_sfs_fs->fp, ",%d,%ld,%e,%e", pair->n_em_iter, pair->snSites, pair->d, args->tole);

	if (args->doDist == 1 && args->do_square_distance == 1)
	{
		fprintf(out_sfs_fs->fp, ",%f,%f", (double)(1.0 - (double)MATH::EST::Sij(pair->SFS)), SQUARE((double)(1.0 - (double)MATH::EST::Sij(pair->SFS))));
	}
	else
	{
		exit(1);
	}
	fprintf(out_sfs_fs->fp, "\n");
}

/// @brief print_SFS_GT print SFS_GT3
/// @param TYPE type of analysis
/// @param out_sfs_fs output file
/// @param args pointer to argStruct
/// @param SFS_GT3 matrix of 3 GT SFS for pair (int **SFS_GT3[pidx])
/// @param snSites (shared) number of sites
void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char *sample1, const char *sample2)
{

	fprintf(out_sfs_fs->fp, "%s,%s,%s,", TYPE, sample1, sample2);
	fprintf(out_sfs_fs->fp, "%d,%d,%d,%d,%d,%d,%d,%d,%d",
			SFS_GT3[0], SFS_GT3[1], SFS_GT3[2],
			SFS_GT3[3], SFS_GT3[4], SFS_GT3[5],
			SFS_GT3[6], SFS_GT3[7], SFS_GT3[8]);

	fprintf(out_sfs_fs->fp, ",%s,%d,%s,%s", TYPE, snSites, TYPE, TYPE);

	if (args->doDist == 1 && args->do_square_distance == 1)
	{
		fprintf(out_sfs_fs->fp, ",%f,%f", (double)(1.0 - (double)MATH::EST::Sij(SFS_GT3, snSites)), SQUARE((double)(1.0 - (double)MATH::EST::Sij(SFS_GT3, snSites))));
		// fprintf(out_sfs_fs->fp, ",%f", SQUARE((double)(1.0 - (double)MATH::EST::Sij(SFS_GT3, snSites))));
	}
	else
	{
		exit(1);
	}
	fprintf(out_sfs_fs->fp, "\n");
}

void DATA::blobStruct_destroy(DATA::blobStruct *c){

	for (size_t i = 0; i < (size_t) c->nContigs; i++)
	{
		FREE(c->contigBlockStartPtrs[i]);
		FREE(c->contigNames[i]);
		FREE(c->contigBlockStarts[i]);
	}
	FREE(c->contigBlockStarts);
	FREE(c->contigNames);
	FREE(c->contigLengths);
	FREE(c->contigNBlocks);
	FREE(c->contigBlockStartPtrs);

	delete c;
}

DATA::blobStruct *DATA::blobStruct_init(const int nContigs, const int blockSize, bcf_hdr_t *hdr){

	DATA::blobStruct *c= new DATA::blobStruct();

	c->nContigs = (size_t)nContigs;
	c->contigNames = (char **)malloc(nContigs * sizeof(char *));
	c->contigLengths = (int *)malloc(nContigs * sizeof(int));

	c->contigBlockStarts = (int **)malloc(nContigs * sizeof(int *));

	c->contigBlockStartPtrs = (double ***)malloc(nContigs * sizeof(double **));
	c->contigNBlocks = (int *)malloc(nContigs * sizeof(int));

	for (size_t i = 0; i < c->nContigs; i++)
	{
		c->contigNames[i] = NULL;
		c->contigBlockStarts[i] = NULL;
		c->contigBlockStartPtrs[i] = NULL;
	}

	for (size_t ci = 0; ci < c->nContigs; ci++)
	{

		const int contigSize = hdr->id[BCF_DT_CTG][ci].val->info[0];
		c->contigLengths[ci] = contigSize;
		fprintf(stderr, "\nContig %ld length:%d\n", ci, contigSize);
		int nBlocks = 0;

		if (blockSize < contigSize)
		{
			nBlocks = (contigSize / blockSize) + 1;
		}
		else
		{
			nBlocks = 1;
			fprintf(stderr, "\nContig %ld is smaller than block size, setting block size to contig size (%d)\n", ci, contigSize);
		}

		// allocate memory for contigBlockStarts
		c->contigBlockStarts[ci] = (int *)malloc(nBlocks * sizeof(int));
		c->contigBlockStartPtrs[ci] = (double **)malloc(nBlocks * sizeof(double *));
		c->contigNBlocks[ci] = nBlocks;

		// fprintf(stderr, "\nContig %d length:%d nBlocks: %d\n", ci, contigSize, nBlocks);
		for (int bi = 0; bi < nBlocks; bi++)
		{
			int blockStart = bi * blockSize;
			c->contigBlockStarts[ci][bi] = blockStart;
			fprintf(stderr, "\nContig %ld block %d starts at %d\n", ci, bi, c->contigBlockStarts[ci][bi]);
		}
	}

	return c;
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



