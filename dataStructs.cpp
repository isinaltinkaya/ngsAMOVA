#include "io.h"
#include "dataStructs.h"

void distanceMatrixStruct::print(IO::outputStruct *out_dm_fs)
{
	fprintf(stderr, "[INFO]\t-> Writing distance matrix to %s.\n", out_dm_fs->fn);

	char buf_values[8192];
	char *bufptr = buf_values;

	for (int px = 0; px < nIndCmb; px++)
	{
		if (px != 0 && px != nIndCmb - 1)
		{
			bufptr += sprintf(bufptr, ",%.*f", (int)DBL_MAXDIG10, M[px]);
		}
		else if (px == 0)
		{
			bufptr += sprintf(bufptr, "%.*f", (int)DBL_MAXDIG10, M[px]);
		}
		else if (px == nIndCmb - 1)
		{
			sprintf(bufptr, ",%.*f\n", (int)DBL_MAXDIG10, M[px]);
		}
		else
		{
			ASSERT(0 == 1);
		}
	}

	switch (out_dm_fs->fc)
	{
	case OUTFC::NONE:
		fprintf(out_dm_fs->fp, "%s", buf_values);
		break;
	case OUTFC::GZ:
		// gzwrite(out_dm_fs->gzfp, M, nIndCmb * sizeof(double));
		gzprintf(out_dm_fs->gzfp, "%s", buf_values);
		break;
	default:
		fprintf(stderr, "[ERROR]\t-> Unknown output file format.\n");
		exit(1);
	}
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
		int sidx = 0;

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

			// ----------------------- INPUT: vcfd ----------------------- //
			// vcf input has pars->nInd and pars->nIndCmb already set from vcfd reading
			if (pars->in_ft == IN_VCF)
			{

				// if input file type is vcf, find the index of the individual in the bcf header
				for (sidx = 0; sidx < pars->nInd + 1; sidx++)
				{

					if (sidx == pars->nInd)
					{
						fprintf(stderr, "\n\n======\n[ERROR] Sample %s not found in the metadata file. \n\n", tok);
						exit(1);
					}
					else if (strcmp(tok, sampleSt->sampleNames[sidx]) == 0)
					{
						pars->vprint(2, "Found sample (vcf_index=%d,id=%s) in the metadata file.", sidx, tok);
						
						// break when found, so that sidx is now the index of the sample in the vcf file
						break;
					}
				}

				// loop through the hierarchical levels in the line
				for (int lvl_i = 0; lvl_i < mS->nLevels; lvl_i++)
				{

					// index of the strata in the current hierarchical level
					int lvl_strata_i = 0;

					// first column after individual is the highest hierarchical level
					tok = strtok(NULL, METADATA_DELIMS);

					// if we are reading the first line
					if (mS->hierArr[lvl_i] == NULL)
					{
						// we are in the first line, initialize the hierStruct for this level

						// TODO
						//  mS->addHierStruct(lvl_i, tok);
						mS->hierArr[lvl_i] = new hierStruct(tok);
						mS->hierArr[lvl_i]->nIndPerStrata[0]++;

						// EOL end of the first line
						if (lvl_i == mS->nLevels - 1)
						{
							// EOL
							// lowest level (excluding individual) in the first line
							// no check because this is the first line of data we read
							ASSERT(lvl_strata_i == 0);
							key = mS->setKeyDigitAtLevel(key, lvl_i, lvl_strata_i);
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

				// sampleSt->addSampleName(sidx, tok);
				sampleSt->addSample(sidx, tok);
				// loop through the hierarchical levels in the line
				for (int lvl_i = 0; lvl_i < mS->nLevels; lvl_i++)
				{

					// index of the strata in the current hierarchical level
					int lvl_strata_i = 0;

					// first column after individual is the highest hierarchical level
					tok = strtok(NULL, METADATA_DELIMS);

					// if we are reading the first line
					if (mS->hierArr[lvl_i] == NULL)
					{
						// we are in the first line, initialize the hierStruct for this level

						// TODO
						//  mS->addHierStruct(lvl_i, tok);
						mS->hierArr[lvl_i] = new hierStruct(tok);
						mS->hierArr[lvl_i]->nIndPerStrata[0]++;

						// EOL end of the first line
						if (lvl_i == mS->nLevels - 1)
						{
							// EOL
							// lowest level (excluding individual) in the first line
							// no check because this is the first line of data we read
							ASSERT(lvl_strata_i == 0);
							key = mS->setKeyDigitAtLevel(key, lvl_i, lvl_strata_i);
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

		if (pars->in_ft == IN_DM)
		{
			pars->nInd = mS->nIndMetadata;
			pars->nIndCmb = nChoose2[pars->nInd];
			pars->init_LUTs();
			set_lut_indsToIdx_2way(pars->nInd, pars->nIndCmb, pars->lut_indsToIdx, pars->lut_idxToInds);
		}
		return (mS);
	}
	else
	{
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
	// fprintf(stderr,"\t--in\t\t\t: input vcfd/BCF filed\n");

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
