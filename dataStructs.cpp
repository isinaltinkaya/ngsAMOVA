#include "io.h"
#include "dataStructs.h"

// TODO is this the proper way to do this?
/// @brief calculate buffer size for CSV output
/// @param n_vals number of values
/// @param max_digits maximum number of digits
/// @return buffer size
const int calculateBufferSizeCsv(const int n_vals, const int max_digits)
{
	// max_digits + 1 = 1 comma
	// result + 2 = 1 newline at the end + 1 null terminator
	return ((n_vals * (max_digits + 1)) + 2);
}

void distanceMatrixStruct::print(IO::outputStruct *out_dm_fs)
{
	fprintf(stderr, "[INFO]\t-> Writing distance matrix to %s.\n", out_dm_fs->fn);
	kstring_t *kbuf = kbuf_init();

	for (int px = 0; px < nIndCmb; px++)
	{
		if (px != 0 && px != nIndCmb - 1)
		{
			// TODO check snprintf maxsize
			//  bufptr += sprintf(bufptr, ",%.*f", (int)DBL_MAXDIG10, M[px]);
			//  bufptr += snprintf(bufptr, max_buf_size, ",%.*f", (int)DBL_MAXDIG10, M[px]);
			//  bufptr += snprintf(bufptr, max_digits+1, ",%.*f", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, ",");
			ksprintf(kbuf, "%.*f", (int)DBL_MAXDIG10, M[px]);
		}
		else if (px == 0)
		{
			// bufptr += sprintf(bufptr, "%.*f", (int)DBL_MAXDIG10, M[px]);
			// bufptr += snprintf(bufptr, max_buf_size, "%.*f", (int)DBL_MAXDIG10, M[px]);
			// bufptr += snprintf(bufptr, max_digits, "%.*f", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, "%.*f", (int)DBL_MAXDIG10, M[px]);
		}
		else if (px == nIndCmb - 1)
		{
			// sprintf(bufptr, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
			// snprintf(bufptr, max_buf_size, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
			// snprintf(bufptr, max_digits+2, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, ",");
			ksprintf(kbuf, "%.*f", (int)DBL_MAXDIG10, M[px]);
			ksprintf(kbuf, "\n");
		}
		else
		{
			NEVER();
		}
	}
	out_dm_fs->write(kbuf);
	kbuf_destroy(kbuf);
}

/// @brief read distance matrix file
/// @param in_dm_fp input distance matrix file
/// @param pars paramStruct parameters
/// @return distance matrix double*
distanceMatrixStruct *distanceMatrixStruct_read_csv(paramStruct *pars, argStruct *args, metadataStruct *metadataSt)
{

	int dm_vals_size = 1225;
	double *dm_vals = (double *)malloc((dm_vals_size) * sizeof(double));

	int n_vals = 0;

	if (IO::isGzFile(args->in_dm_fn) == 1)
	{

		size_t buf_size = FGETS_BUF_SIZE;
		size_t *buf_size_ptr = &buf_size;
		char *line = (char *)malloc(buf_size);
		char **line_ptr = &line;
		IO::readGzFile::readToBuffer(args->in_dm_fn, line_ptr, buf_size_ptr);

		ASSERT(line != NULL);
		char *tok = strtok(line, ",\n");
		while (tok != NULL)
		{
			if (n_vals > dm_vals_size)
			{
				dm_vals_size = dm_vals_size * 2;
				dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
			}
			dm_vals[n_vals] = atof(tok);
			n_vals++;
			tok = strtok(NULL, ",\n");
		}

		FREE(line);
		FREE(tok);

		// while (gzgets(fp, line, buf_size) != NULL)
		// {
		// 	// check if the line was fully read
		// 	size_t line_len = strlen(line);
		// 	if (line[line_len - 1] == '\n')
		// 	{
		// 		// line was fully read
		// 		break;
		// 	}
		// 	else
		// 	{
		// 		fprintf(stderr, "\t-> Line was not fully read, increasing buffer size\n");
		// 		// line was not fully read
		// 		buf_size *= 2;

		// 		char *new_line = new char[buf_size];
		// 		new_line = (char *)realloc(line, buf_size);
		// 		ASSERT(new_line != NULL);
		// 		line = new_line;
		// 	}
		// }
		// gzclose(fp);
	}
	else
	{

		int buf_size = IO::readFile::getBufferSize(args->in_dm_fn);
		char *line = (char *)malloc(buf_size);
		ASSERT(line != NULL);

		char dm_buf[buf_size];

		FILE *in_dm_fp = fopen(args->in_dm_fn, "r");
		while (fgets(dm_buf, buf_size, in_dm_fp))
		{
			char *tok = strtok(dm_buf, ",\n");
			while (tok != NULL)
			{
				if (n_vals > dm_vals_size)
				{
					dm_vals_size = dm_vals_size * 2;
					dm_vals = (double *)realloc(dm_vals, (dm_vals_size) * sizeof(double));
				}
				dm_vals[n_vals] = atof(tok);
				n_vals++;
				tok = strtok(NULL, ",\n");
			}
		}

		FREE(line);
		FCLOSE(in_dm_fp);
	}

	fprintf(stderr, "[INFO]\t-> Number of values in distance matrix: %d. (i.e. number of unique individual pairs)\n", n_vals);
	fprintf(stderr, "[INFO]\t-> Number of individuals based on the number of individuals in the distance matrix: %d.\n", find_n_given_nC2(n_vals));

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

	FREE(dm_vals);
	return dMS;
}

// if formula is not defined, hierarchical levels must defined in this order:
// individual, highest_level, ..., lowest_level
// e.g.
// individual, region, population, subpopulation
// which translates to the formula:
// individual ~ region / population / subpopulation
metadataStruct *metadataStruct_get(FILE *in_mtd_fp, sampleStruct *sampleSt,
								   formulaStruct *FORMULA, int has_colnames,
								   paramStruct *pars)
{

	char *mtd_buf = new char[FGETS_BUF_SIZE];

	size_t nLevels = 0;

	// go to beginning of file
	ASSERT(fseek(in_mtd_fp, 0, SEEK_SET) == 0);

	if (FORMULA == NULL)
	{

		if (has_colnames != 1)
		{
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
			IO::vprint(1, "\n[INFO]\t-> Found hierarchical level: %s\n", hdrtok);
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
		int sidx = 0;

		// loop through the remaining lines; one line per individual
		while (fgets(mtd_buf, FGETS_BUF_SIZE, in_mtd_fp))
		{

			mS->nInd++;

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
			// vcf input has pars->nInd and pars->nIndCmb already set from vcfd reading
			if (pars->in_ft == IN_VCF)
			{

				// if input file type is vcf, find the index of the individual in the bcf header
				// TODO check nInd+1 why
				for (sidx = 0; sidx < pars->nInd + 1; sidx++)
				{

					if (sidx == pars->nInd)
					{
						fprintf(stderr, "\n\n======\n[ERROR] Sample %s not found in the metadata file. \n\n", tok);
						exit(1);
					}
					else if (strcmp(tok, sampleSt->sampleNames[sidx]) == 0)
					{
						IO::vprint(2, "Found sample (vcf_index=%d,id=%s) in the metadata file.", sidx, tok);

						// break when found, so that sidx is now the index of the sample in the vcf file
						break;
					}
				}

				// index of the strata in the current hierarchical level
				int strata_idx_i = 0;
				// loop through the hierarchical levels in the line
				for (int lvl_i = 0; lvl_i < mS->nLevels; lvl_i++)
				{

					// first column after individual name is the highest hierarchical level
					tok = strtok(NULL, METADATA_DELIMS);

					// if we are reading the first line
					if (mS->hierArr[lvl_i] == NULL)
					{
						// we are in the first line, initialize the hierStruct for this level

						mS->hierArr[lvl_i] = new hierStruct(lvl_i);
						strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						ASSERT(strata_idx_i == 0); // every index in the first line should be 0

						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
						// add the strata to its parent's subStrataIdx
						mS->addSubStrataIdx(key, lvl_i);

						// EOL end of the first line
						if (lvl_i == mS->nLevels - 1)
						{
							// EOL
							// lowest level (excluding individual) in the first line
							mS->ind2stratakey[sidx] = key;
							continue;
						}
					}

					// EOL we are at the lowest level column (of a line that is not the first line)
					// associate the individual with the stratakey
					else if (lvl_i == mS->nLevels - 1)
					{
						strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
						// add the strata to its parent's subStrataIdx
						mS->addSubStrataIdx(key, lvl_i);
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

						strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
						// add the strata to its parent's subStrataIdx
						mS->addSubStrataIdx(key, lvl_i);
					}
				}
			}
			// ----------------------- INPUT: DM ----------------------- //
			else if (pars->in_ft == IN_DM)
			{

				// index of the strata in the current hierarchical level
				int strata_idx_i = 0;

				sampleSt->addSample(sidx, tok);
				// loop through the hierarchical levels in the line
				for (int lvl_i = 0; lvl_i < mS->nLevels; lvl_i++)
				{

					// first column after individual is the highest hierarchical level
					tok = strtok(NULL, METADATA_DELIMS);

					// if we are reading the first line
					if (mS->hierArr[lvl_i] == NULL)
					{
						// we are in the first line, initialize the hierStruct for this level
						mS->hierArr[lvl_i] = new hierStruct(lvl_i);
						strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						ASSERT(strata_idx_i == 0);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
						// add the strata to its parent's subStrataIdx
						mS->addSubStrataIdx(key, lvl_i);

						// EOL end of the first line
						if (lvl_i == mS->nLevels - 1)
						{
							// EOL
							// lowest level (excluding individual) in the first line
							mS->ind2stratakey[sidx] = key;
							sidx++; // TODO not in vcf
							continue;
							// }else if(lvl_i != 0)
							// {
						}
					}

					// EOL
					// at the lowest level column; end of the individual's line
					else if (lvl_i == mS->nLevels - 1)
					{
						// associate the individual with the stratakey

						strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
						// add the strata to its parent's subStrataIdx
						mS->addSubStrataIdx(key, lvl_i);
						mS->ind2stratakey[sidx] = key;

						sidx++; // TODO
					}
					else
					{

						// at a regular column, nothing special about it
						// (not the first line; not the first column; not the last column)
						//
						// e.g. ind1,x1,y1,z1
						// 		ind2,[x1],y2,z3
						//   		we are at 'x1' of ind2, already allocated in previous individual

						strata_idx_i = mS->hierArr[lvl_i]->getStrataIndex(tok);
						key = mS->setKeyDigitAtLevel(key, lvl_i, strata_idx_i);
						// add the strata to its parent's subStrataIdx
						mS->addSubStrataIdx(key, lvl_i);
					}
				}
			}
			else
			{
				fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
				ASSERT(0 == 1);
			}
		}
		// mS->resize();

		// TODO printtable
		//  mS->print(stderr);

		// TODO improve with bitarray LUT, 1 if pair belongs to given strata 0 otherwise
		//  pairToAssoc[pair_idx][hier_level][strata_idx] = 1
		for (size_t i = 0; i < (size_t)mS->nLevels + 1; i++)
		{
			delete[] levelNames[i];
		}
		delete[] levelNames;
		delete[] mtd_buf;

		pars->nInd = mS->nInd;
		pars->nIndCmb = NC2_LUT[pars->nInd];
		mS->nIndCmb = pars->nIndCmb;
		if (pars->in_ft == IN_DM)
		{
			pars->init_LUTs();
			set_lut_indsToIdx_2way(pars->nInd, pars->nIndCmb, pars->lut_indsToIdx, pars->lut_idxToInds);
		}

		ASSERT(mS->nInd > 0);
		mS->pairToAssoc_create();
		for (int l = 0; l < mS->nLevels; l++)
		{
			for (int gi = 0; gi < mS->hierArr[l]->nStrata; gi++)
			{
				for (int i1 = 0; i1 < pars->nInd - 1; i1++)
				{
					for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
					{
						if (mS->pairInStrataAtLevel(i1, i2, l, gi) == 1)
						{
							mS->pairToAssoc[pars->lut_indsToIdx[i1][i2]][l][gi] = 1;
						}
					}
				}
			}
		}

		return (mS);
	}
	else
	{ // formula is defined, i.e. not NULL
		fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
		ASSERT(0 == 1);
	}

	NEVER();
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
	// fprintf(fp, "Version: %s (using htslib %s)\n\n", program_version(), hts_version());
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
			" -i\n"
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
	exit(0);
}

/// @brief formulaStruct_get initialize the formulaStruct
/// @param formula formula string
/// @return pointer to formulaStruct
/// @example formula = 'Samples ~ Continents/Regions/Populations'
formulaStruct *formulaStruct_get(const char *formula)
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
// void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, pairStruct *pair, argStruct *args, const char *sample1, const char *sample2)
// {

// 	fprintf(out_sfs_fs->fp, "%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f",
// 			TYPE,
// 			sample1,
// 			sample2,
// 			pair->snSites * pair->SFS[0], pair->snSites * pair->SFS[1], pair->snSites * pair->SFS[2],
// 			pair->snSites * pair->SFS[3], pair->snSites * pair->SFS[4], pair->snSites * pair->SFS[5],
// 			pair->snSites * pair->SFS[6], pair->snSites * pair->SFS[7], pair->snSites * pair->SFS[8]);

// 	fprintf(out_sfs_fs->fp, ",%d,%ld,%e,%e", pair->n_em_iter, pair->snSites, pair->d, args->tole);

// 	if (args->doDist == 1 && args->do_square_distance == 1)
// 	{
// 		fprintf(out_sfs_fs->fp, ",%f,%f", (double)(1.0 - (double)MATH::EST::Sij(pair->SFS)), SQUARE((double)(1.0 - (double)MATH::EST::Sij(pair->SFS))));
// 	}
// 	else
// 	{
// 		exit(1);
// 	}
// 	fprintf(out_sfs_fs->fp, "\n");
// }

/// @brief print_SFS_GT print SFS_GT3
/// @param TYPE type of analysis
/// @param out_sfs_fs output file
/// @param args pointer to argStruct
/// @param SFS_GT3 matrix of 3 GT SFS for pair (int **SFS_GT3[pidx])
/// @param snSites (shared) number of sites
// void IO::print::Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char *sample1, const char *sample2)
// {

// 	fprintf(out_sfs_fs->fp, "%s,%s,%s,", TYPE, sample1, sample2);
// 	fprintf(out_sfs_fs->fp, "%d,%d,%d,%d,%d,%d,%d,%d,%d",
// 			SFS_GT3[0], SFS_GT3[1], SFS_GT3[2],
// 			SFS_GT3[3], SFS_GT3[4], SFS_GT3[5],
// 			SFS_GT3[6], SFS_GT3[7], SFS_GT3[8]);

// 	fprintf(out_sfs_fs->fp, ",%s,%d,%s,%s", TYPE, snSites, TYPE, TYPE);

// 	if (args->doDist == 1)
// 	{
// 		if (args->do_square_distance == 1)
// 		{
// 			fprintf(out_sfs_fs->fp, "%f", (double)SQUARE(MATH::EST::Dij(SFS_GT3, snSites)));
// 		}
// 		else
// 		{
// 			fprintf(out_sfs_fs->fp, "%f", (double)MATH::EST::Dij(SFS_GT3, snSites));
// 		}
// 	}
// 	else
// 	{
// 		exit(1);
// 	}
// 	fprintf(out_sfs_fs->fp, "\n");
// }

void blobStruct_destroy(blobStruct *c)
{

	for (size_t i = 0; i < (size_t)c->nContigs; i++)
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

blobStruct *blobStruct_init(const int nContigs, const int blockSize, bcf_hdr_t *hdr)
{

	blobStruct *c = new blobStruct();

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

// /// Warnings
// ///
// /// Warnings[IsRelatedToTypeX][IndexInTypeX] = "Warning string"
// ///
// /// Why do you get this warning = You used [IsRelatedToTypeX][IndexInTypeX]
// /// e.g. Warnings[INFT][IN_VCF] = "You used VCF input file type"
// const char* WARNINGS[][] = {
// 	// INFT input file type
// 	{
// 		// IN_VCF	input file type is VCF
// 		{"Assuming that the VCF file is sorted by position."},
// 		// IN_DM	input file type is distance matrix
// 		{"Assuming that the distance matrix is either an output from this program or a distance matrix with prepared with the same format as the output of this program."},
// 		// IN_JGPD
// 		{""}
// 	}
// }

/// @brief  getDigitIndexAtLevel - get the index of a digit from the right
/// @param nLevels  			total number of levels
/// @param lvl  				0-indexed level
/// @param maxDigitsPerHLevel  	maximum number of digits allowed in key construction per hierarchical level
/// @return  					index of digit (from the right)
/// @details  					Example: nLevels = 3, lvl = 1, maxDigitsPerHLevel = 2
///								[1]	[][] [][] [][]
///								 ^ 	  ^    ^    ^
///								 ^	  ^	   ^	|______ the subkey for lvl 2 starts from here (lowest level, e.g. subpopulation)
/// 							 ^	  ^	   |______ the subkey for lvl 1 starts from here (second highest level, e.g. population)
///								 ^    |______ the subkey for lvl 0 starts from here (highest level, e.g. region)
///								 |______ the key holder starts from here (holding the rest together)
///
int getDigitIndexAtLevel(const int nLevels, const int lvl, const int maxDigitsPerHLevel)
{

	// convert to 0-indexed
	int n = nLevels - 1;
	return (n - lvl) * maxDigitsPerHLevel;
}

/// @brief   [overload] getDigitIndexAtLevel - get the index of a digit from the right
/// @details function overload for default MAXDIG_PER_HLEVEL
int getDigitIndexAtLevel(const int nLevels, const int lvl)
{
	int n = nLevels - 1;
	return (n - lvl) * MAXDIG_PER_HLEVEL;
}

/// @brief calculateKeyAtLevel - calculate the strata key value at a given strata index in a given level
/// @param lvl			- hierarchical level
/// @param strata_idx 	- strata index
/// @return integer value of the strata key
/// @example lvl=1, strata_idx=2 -> 1[02]00
int calculateKeyAtLevel(const int lvl, const int strata_idx)
{
	return (int)(pow(10, MAXDIG_PER_HLEVEL * (lvl + 1)) + strata_idx);
	// POW10_LUT[]
}

/// @brief estimate_dxy - estimate dxy statistic for a pair of individuals
/// @param idx1			- index of the first group at specified level
/// @param idx2			- index of the second group at specified level
/// @param lvl			- hierarchical level the group indices refer to
/// @param dMS          - distance matrix struct
/// @param mtd          - metadata struct
/// @param pars         - parameter struct
/// @return dxy         - dxy statistic
double estimate_dxy(const int idx1, const int idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars)
{

	double dxy = 0.0;

	// TODO below
	//  // extract the strata id at level from key1 and key2
	//  const int strata1 = (key1 >> (lvl*8)) & 0xFF;
	//  const int strata2 = (key2 >> (lvl*8)) & 0xFF;

	// lvl is 0-indexed and nLevels is count
	if (lvl >= mtd->nLevels)
	{
		fprintf(stderr, "\n[ERROR][estimate_dxy] The level specified (%d) is greater than the number of levels (%d)\n", lvl + 1, mtd->nLevels);
		exit(1);
	}

	// get number of individuals belonging to strata1 (extracted from key1 based on lvl)
	const int nInd1 = mtd->hierArr[lvl]->nIndPerStrata[idx1];
	const int nInd2 = mtd->hierArr[lvl]->nIndPerStrata[idx2];

	IO::vprint(1, "Running estimate_dxy for strata %d (nInd:%d) and strata %d (nInd:%d) at level %d", idx1, nInd1, idx2, nInd2, lvl);

	double nxny = (double)(nInd1 * nInd2);

	if (idx1 == idx2)
	{
		fprintf(stderr, "\n[ERROR][estimate_dxy] idx1:%d is equal to idx2:%d\n", idx1, idx2);
		exit(1);
	}

	// locate the individual pairs in the distance matrix where one individual is from group 1 and the other is from group 2
	for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
	{

		// use i1s belonging to idx1
		if (mtd->indInStrata(i1, lvl, idx1) != 1)
			continue;

		for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
		{

			// use i2s belonging to idx2
			if (mtd->indInStrata(i2, lvl, idx2) != 1)
				continue;

			IO::vprint(2, "Running estimate_dxy for i1:%d from group with index %d and i2:%d from group with index %d at hierarchical level %d\n", i1, idx1, i2, idx2, lvl);

			// IO::vprint(3,"dxy was %f and adding %f",dxy,dMS->M[pars->lut_indsToIdx[i1][i2]]);
			dxy += dMS->M[pars->lut_indsToIdx[i1][i2]];
		}
	}
	dxy = dxy / nxny;
	return dxy;
}