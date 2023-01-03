#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include <limits>
#include <cstddef>
#include <time.h>

#include <cstdlib>
#include <cstdio>

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <math.h>

#ifndef MAX_FORMULA_TOKENS
#define MAX_FORMULA_TOKENS 10
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/*
 * Macro:[DBL_MAXDIG10]
 * Defines the maximum number of decimal digits that can be represented by a
 * double-precision floating-point number on a 64-bit system
 *
 * Calculates the max number of decimal digits that a double-precision
 * floating-point number can have on a 64-bit system using the value of DBL_MANT_DIG
 * The UL suffix ensures that our constant values are unsigned long
 * to prevent any overflow issues.
 *
 */
#define DBL_MAXDIG10 (2 + (DBL_MANT_DIG * 30103UL) / 100000UL)

/*
 * Macro:[AT]
 * Injects the file and line info as string
 */
#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)

/*
 * Macro:[ASSERT]
 * shortcut to evaluate an expression, works the same way as the C-macro assert
 */
#define ASSERT(expr)                                                                             \
	if (!(expr))                                                                                 \
	{                                                                                            \
		fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); \
		exit(1);                                                                                 \
	}

/*
 * Macro:[FREE]
 * shortcut to free memory and set pointer to NULL
 * and catch double free
 */
#define FREE(expr)                                                                                             \
	if (expr)                                                                                                  \
	{                                                                                                          \
		free(expr);                                                                                            \
		expr = NULL;                                                                                           \
	}                                                                                                          \
		// fprintf(stderr, "\n\n*******\n[FREEING NULL MEMORY](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); 


extern const int get_3x3_idx[3][3];

const double NEG_INF = -std::numeric_limits<double>::infinity();

using size_t = decltype(sizeof(int));

/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_fn		pointer to input file name
 * @field *in_mtd_fn	pointer to input Metadata file name
 * @field *out_fn		pointer to output file prefix [angsdput]
 * @field seed			random seed
 *
 * @field isSim			input is vcfgl simulation output
 * 							anc=ref and der=alt[0]
 *
 * @field isTest		[DEV] run for testing purposes
 *
 * @field minInd		[-1 = not set]
 * 						minimum number of individuals needed
 * 						for site to be included in analyses
 *
 * 						if minInd not set; set minInd=2
 * 						== no filter, include all sites that exist in pair
 *
 * 						if minInd==nInd; set minInd=0
 * 						== site should be nonmissing for all individuals
 * 						to be included
 *
 * @field blockSize		[0 = not set]
 *						Block size to be used in block bootstrapping
 *
 *
 * @field formula		[not set]
 * 						the formula of the AMOVA model to be fitted
 *
 * 						of the form:
 *  						{LEFT_HAND_SIDE} ~ {RIGHT_HAND_SIDE}
 *
 *	 						LEFT_HAND_SIDE: Column defining the distance matrix
 *											i.e. Samples
 *
 *	 						RIGHT_HAND_SIDE: Column(s) defining the hierarchical stratification levels
 *											i.e. Populations, Regions, etc
 *
 *
 * 						e.g. 1 level
 * 							Samples ~ Populations
 * 							{SamplesColumn} ~ {HierLvl1_StrataColumn}
 *
 * 						e.g. 2 levels
 * 							Samples ~ Regions/Populations
 * 							{SamplesColumn} ~ {HierLvl1_StrataColumn} / {HierLvl2_StrataColumn}
 *
 * 						e.g. 3 levels
 * 							Samples ~ Continents/Regions/Populations
 * 							{SamplesColumn} ~ {HierLvl1_StrataColumn} / {HierLvl2_StrataColumn} / {HierLvl3_StrataColumn}
 *
 * 						NOTE: All column names must exist in Metadata file
 * 						e.g.
 *
 * 						Samples,Populations,Regions
 * 						ind1,Pop1,Reg1
 * 						ind2,Pop1,Reg1
 * 						ind3,Pop2,Reg2
 * 						ind4,Pop2,Reg2
 *
 *
 * @field keyCols		defines the index(es) (1-based) of the
 * 						column(s) in Metadata file 'in_mtd_fn'
 * 						to be used to define hierarchical stratification levels
 *
 * @field doAMOVA		[0]
 * 						1 use 10 genotype likelihoods (GL)
 * 						2 use genotypes (GT) (NOTE: Only for benchmark purposes)
 *
 * @field doDist		[0] use Sij similarity index
 * 						[1] use Dij (1-Sij) dissimilarity index
 * 						[2] use Fij F statistic [DEPRECATED]
 *
 * @field sqDist		[0] use absolute value of distance measure (|dist_ij|)
 * 						[1] use squared distance measure (dist_ij^2)
 *
 * @field doInd			do ind pairs
 * @field ind1			ind1 id
 * @field ind2			ind2 id
 *
 *
 * @field doTest		test for convergence
 *
 *
 * @field mThreads		maximum number of threads defined by user
 * @field mEmIter		maximum number of iterations allowed for em
 *
 */
typedef struct
{

	char *in_fn;
	char *in_sfs_fn;
	char *in_mtd_fn;
	char *out_fn;

	char *formula;
	int *keyCols;
	int hasColNames;

	int blockSize;
	int doAMOVA;
	int printMatrix;
	int isSim;
	int isTest;
	int doDist;
	int sqDist;
	int minInd;

	int seed;

	double tole;

	int doInd;
	int ind1;
	int ind2;

	int doTest;

	int mThreads;
	int mEmIter;

	int gl2gt;

} argStruct;

argStruct *argStruct_init();
argStruct *argStruct_get(int argc, char **argv);
// void *argStruct_destroy(argStruct *arg);

/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites				number of sites
 * @field nInd					number of individuals
 * @field pos					position
 *
 * @field LUT_indPair_idx		lookup table for mapping two individuals
 * 								to their pair index
 * @field n_ind_cmb				number of unique pairwise individual combinations
 *
 * @field major					major allele
 * @field minor					minor allele
 * @field ref					reference allele
 * @field anc					ancestral allele
 * @field der					derived allele
 */
typedef struct
{

	size_t nSites;
	size_t totSites;

	int nInd;

	int *pos;

	int **LUT_indPair_idx;
	int n_ind_cmb;

	char *major;
	char *minor;
	char *ref;
	char *anc;
	char *der;

	char *DATETIME;

	void print()
	{
		// print lookup table
		for (int i1 = 0; i1 < nInd - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < nInd; i2++)
			{
				fprintf(stderr, "\n%i %i %i\n", LUT_indPair_idx[i1][i2], i1, i2);
			}
		}
	}

} paramStruct;

paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

char *get_time();

void usage(FILE *fp);

namespace DATA
{
	typedef struct formulaStruct
	{
		int nTokens;
		int nTokensFound = 0;
		const char *formula = NULL;
		char **formulaTokens;
		size_t *formulaTokenIdx;

		void print(FILE *fp)
		{
			fprintf(fp, "\nFormula: %s", formula);
			fprintf(fp, "\nTokens: %i", nTokens);
			for (int i = 0; i < nTokens; i++)
			{
				// fprintf(fp, "\n\t%i\t%s\t%i", i, formulaTokens[i], formulaTokenIdx[i]);
				fprintf(fp, "\n\t%i\t%s\n", i, formulaTokens[i]);
			}
		}

	} formulaStruct;
	formulaStruct *formulaStruct_get(const char *formula);

	typedef struct contigsStruct
	{

		// Number of contigs
		size_t nContigs;

		// array of contig indices
		int *contigIdx;

		// 2d array of contig names
		// [Contig][Name]
		char **contigNames;

		// Array of number of blocks per contig
		int *contigNBlocks;

		// Array of contig lengths
		int *contigLengths;

		// 2D array of contig block starts
		// [Contig][BlockStart]
		// size_t *contigBlockStarts;
		int **contigBlockStarts;

		// [Contig][Block] int* to block start
		// 2D array of pointers to actual contig block starts
		double ***contigBlockStartPtrs;

		// void contig_set_ContigBlockStarts(int ci, int nBlocks){
		// 	contigBlockStarts[ci] = (int *)realloc(contigBlockStarts[ci], nBlocks * sizeof(int));
		// }

		contigsStruct(const int nContigs_)
		{

			nContigs = (size_t)nContigs_;
			contigNames = (char **)malloc(nContigs * sizeof(char *));
			contigLengths = (int *)malloc(nContigs * sizeof(int));

			contigBlockStarts = (int **)malloc(nContigs * sizeof(int *));

			contigBlockStartPtrs = (double ***)malloc(nContigs * sizeof(double **));
			contigNBlocks = (int *)malloc(nContigs * sizeof(int));

			for (size_t i = 0; i < nContigs; i++)
			{
				contigNames[i] = NULL;
				contigBlockStarts[i] = NULL;
				contigBlockStartPtrs[i] = NULL;
			}
		}

		~contigsStruct()
		{
			for (size_t i = 0; i < nContigs; i++)
			{

				FREE(contigBlockStartPtrs[i]);

				FREE(contigNames[i]);

				FREE(contigBlockStarts[i]);
			}
			FREE(contigBlockStarts);
			FREE(contigNames);
			FREE(contigLengths);
			FREE(contigNBlocks);

			FREE(contigBlockStartPtrs);
			contigBlockStartPtrs = NULL;
		}

	} contigsStruct;
	contigsStruct *contigsStruct_init(const int n_contigs);

	typedef struct pairStruct
	{

		int i1;
		int i2;
		int idx;

		size_t snSites = 0;

		int *sSites;
		// TODO init better
		size_t _sSites = 100;

		// int nDim;

		double d;
		int n_em_iter;

		// TODO should not always create below
		//  double SFS[3][3]={{NEG_INF,NEG_INF,NEG_INF},
		//  {NEG_INF,NEG_INF,NEG_INF},
		//  {NEG_INF,NEG_INF,NEG_INF}};

		double *SFS = NULL;

		pairStruct(int ind1, int ind2, int pair_idx)
		{
			i1 = ind1;
			i2 = ind2;
			idx = pair_idx;
			d = 0.0;
			n_em_iter = 0;
			sSites = (int *)malloc(_sSites * sizeof(int));
			for (size_t i = 0; i < _sSites; i++)
			{
				sSites[i] = -1;
			}

			SFS = (double *)malloc(9 * sizeof(double));
			for (int i = 0; i < 9; i++)
			{
				SFS[i] = NEG_INF;
			}
		}
		~pairStruct()
		{
			FREE(sSites);

			FREE(SFS);
		}

		void print(FILE *fp, bcf_hdr_t *hdr)
		{

			fprintf(fp, "%d,%d,%d,%s,%s", i1, i2, idx, hdr->samples[i1], hdr->samples[i2]);
		}

		void print(FILE *fp)
		{

			fprintf(fp, "%d,%d,%d", i1, i2, idx);
		}

	} pairStruct;

	typedef struct samplesStruct
	{

		int *sampleArr = NULL;
		char **bcfSamples = NULL;
		char **sampleNames = NULL;
		int *sampleOrder = NULL;

		int nSamples = 0;

		samplesStruct(int nSamples_)
		{
			nSamples = nSamples_;
			sampleArr = (int *)malloc(nSamples * sizeof(int));
			bcfSamples = (char **)malloc(nSamples * sizeof(char *));
			sampleNames = (char **)malloc(nSamples * sizeof(char *));
			sampleOrder = (int *)malloc(nSamples * sizeof(int));

			for (int i = 0; i < nSamples; i++)
			{
				sampleArr[i] = 0;
				bcfSamples[i] = NULL;
				sampleNames[i] = NULL;
				sampleOrder[i] = 0;
			}
		}

		~samplesStruct()
		{
			FREE(sampleArr);
			FREE(sampleOrder);

			for (int i = 0; i < nSamples; i++)
			{
				FREE(bcfSamples[i]);

				FREE(sampleNames[i]);
			}
			// TODO
			free(bcfSamples);
			free(sampleNames);
			bcfSamples = NULL;
			sampleNames = NULL;
		}

	} samplesStruct;

	// strata: defines the population factor of the data
	//
	// 		has 1 entry per individual
	// 		corresponding to which strata it belongs to
	//		at the specific hierarchical level
	//
	// 		eg. if there are 2 stratas (hLevel=2), and 3 individuals
	// 		ind,Region,Population
	// 		1,1,1
	// 		2,1,2
	// 		3,2,3
	//
	//		hierStruct[hLevels] = {Region,Population}
	//							= {strataStruct*, strataStruct*}

	//		at hLevel=0, hierStruct[0]
	//			hierStruct[0]->strataStruct->id="Region"
	//			hierStruct[0]->strataStruct->nStrata=2  //num unique values
	//			hierStruct[0]->strataStruct->strataArr = {1,1,2}
	//
	//		at hLevel=1, hierStruct[1]
	//			hierStruct[1]->strataStruct->id="Population"
	//			hierStruct[1]->strataStruct->nStrata=3  //num unique values
	//			hierStruct[1]->strataStruct->strataArr = {1,2,3}
	//
	//
	// strata is defined based on the definition in adegenet
	// https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-strata.pdf
	// typedef struct strataStruct
	// {

	// 	// number of unique values at this hierarchical level
	// 	size_t nStrata;

	// 	char *id;
	// 	size_t _id = 100; // TODO increase if overflow

	// 	strataStruct(size_t nStrata_, char *id_)
	// 	{
	// 		nStrata = nStrata_;

	// 		id = (char *)malloc(_id * sizeof(char));
	// 		strcpy(id, id_);
	// 	}
	// 	~strataStruct()
	// 	{
	// 		free(id);
	// 		id = NULL;
	// 	}

	// } strataStruct;

	// TODO old, new above
	typedef struct strataStruct
	{

		int nInds = 0;
		char *id;

		// TODO Associate hierarchical levels
		//
		// int assoc = 0;

		strataStruct()
		{
			id = NULL;
		}
		~strataStruct()
		{
			free(id);
			id = NULL;
		}

	} strataStruct;

	// Metadata
	// |
	// -- nLevels
	// -- strataArr 1 strataArr

	// strataArr[ind][hLevel][nStrataInLevel]
	// strataArr[ind][hLevel_i]=*strataStruct_i;
	typedef struct metadataStruct
	{

		size_t nLevels = 1;

		// samplesStruct *samples;

		// char **sampleNames = NULL;

		// map hdr->samples to sample names in metadata
		// eg. hdr->samples[0] = MTD->sampleNames[1]
		// 		hdr->samples[1] = MTD->sampleNames[0]
		// 		then access individual 1 using MTD->sampleNames[MTD->sampleOrder[1]]
		//
		// sampleOrder[X] = Y
		// 		X = index in hdr->samples
		// 		Y = index in MTD, the order found in the metadata file
		// int *sampleOrder = NULL;

		int nStrata=0;
		strataStruct *strataArr;

		int nInds_total = 0;

		int* nStrataAtLevel= NULL;
		int nSamples = 0;

		metadataStruct(int nLevels_)
		{
			
			fprintf(stderr, "\n-> Creating metadataStruct");
			nLevels = nLevels_;

			strataArr = new strataStruct[nLevels];
			nStrataAtLevel = new int[nLevels];
			
		};
		~metadataStruct()
		{
			delete[] strataArr;

		}

	} metadataStruct;
	metadataStruct *metadataStruct_get(FILE *in_mtd_fp,
						   int *keyCols, const char *delims, DATA::samplesStruct *SAMPLES,
						   DATA::formulaStruct *FORMULA, int HAS_COLNAMES);
	// //TODO old, new above
	// 	typedef struct metadataStruct
	// 		{

	// 			strataStruct *strataArr;
	// 			size_t _strataArr = 10;

	// 			int nInds_total = 0;

	// 			int nStrata = 0;

	// 			metadataStruct()
	// 			{
	// 				strataArr = new strataStruct[_strataArr];
	// 				fprintf(stderr, "\n->Creating metadataStruct");
	// 			};
	// 			~metadataStruct()
	// 			{
	// 				delete[] strataArr;
	// 				strataArr = NULL;
	// 			}

	// 		} metadataStruct;

};

namespace IO
{

	char *setFileName(const char *a, const char *b);

	FILE *getFile(const char *fname, const char *mode);
	FILE *openFile(const char *a, const char *b);
	FILE *openFile(char *c);

	namespace readFile
	{
		// int Metadata(DATA::metadataStruct *MTD, FILE *in_mtd_fp, int *keyCols, const char *delims, DATA::samplesStruct *SAMPLES, DATA::formulaStruct *FORMULA, int HAS_COLNAMES);
		int SFS(FILE *in_sfs_fp, const char *delims, DATA::samplesStruct *SAMPLES);

		char *getFirstLine(char *fn);
		char *getFirstLine(FILE *fp);
	};

	namespace inspectFile
	{
		int count_nColumns(char *line, const char *delims);
		int count_nRows(char *fn, int HAS_COLNAMES);
		int count_nRows(FILE *fp, int HAS_COLNAMES);
	};

	namespace validateFile
	{
		int Metadata(FILE *in_mtd_fp, int nInds, int *keyCols, DATA::formulaStruct *FORMULA, const char *delims, int HAS_COLNAMES);

	};

	typedef struct outputStruct
	{
		char *fn = NULL;
		FILE *ff = NULL;

		outputStruct(const char *fp, const char *suffix)
		{
			fn = setFileName(fp, suffix);
			ff = openFile(fn);
		}
		~outputStruct()
		{
			fclose(ff);
			free(fn);
			fn = NULL;
		}

	} outputStruct;

	typedef struct outFilesStruct
	{
		outputStruct *out_emtest_fs = NULL;
		outputStruct *out_sfs_fs = NULL;
		outputStruct *out_dm_fs = NULL;
		outputStruct *out_amova_fs = NULL;

		outFilesStruct(argStruct *args)
		{
			if (args->printMatrix == 1)
			{
				out_dm_fs = new outputStruct(args->out_fn, ".distance_matrix.csv");
			}
			if (args->doTest == 1)
			{
				out_emtest_fs = new outputStruct(args->out_fn, ".emtest.csv");
			}

			out_amova_fs = new outputStruct(args->out_fn, ".amova.csv");
			out_sfs_fs = new outputStruct(args->out_fn, ".sfs.csv");
		}

		~outFilesStruct()
		{
			delete out_emtest_fs;
			delete out_dm_fs;
			delete out_sfs_fs;
			delete out_amova_fs;
		}

	} outFilesStruct;

	namespace print
	{

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
		void Array(FILE *fp, double *arr, size_t N, size_t M, char sep);

		// IO::print::Array
		// :overload: int array
		void Array(FILE *fp, int *arr, size_t N, size_t M, char sep);

		/// @brief print sfs output amovaput.sfs.csv
		/// @param TYPE type of analysis
		/// @param out_sfs_fs pointer to output file
		/// @param pair pair of samples
		/// @param pars parameters struct
		/// @param args arguments struct
		/// @param sample1 name of sample 1
		/// @param sample2 name of sample 2
		void Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, DATA::pairStruct *pair, argStruct *args, const char *sample1, const char *sample2);

		void Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);

		void M_PWD(const char *TYPE, IO::outputStruct *out_dm_fs, int n_ind_cmb, double *M_PWD);
	}
}

typedef struct threadStruct
{

	DATA::pairStruct *pair;
	double **lngls;

	FILE *out_sfs_fp;

	// double* M_PWD_GL_PAIR;

	// TODO use them globally?
	double tole;
	int mEmIter;

	size_t nSites;

	threadStruct(DATA::pairStruct *tPair, double **lngl, size_t nSites_t, IO::outputStruct *out_sfs_fs, double toleArg, int mEmIterArg)
	{
		pair = tPair;
		lngls = lngl;
		nSites = nSites_t;
		out_sfs_fp = out_sfs_fs->ff;
		// M_PWD_GL_PAIR=M_PWD_GL_P;
		tole = toleArg;
		// doDist=doDistArg;
		mEmIter = mEmIterArg;
	}

} threadStruct;
void print_SFS_GT(const char *TYPE, IO::outputStruct *out_sfs_fs, paramStruct *pars, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);

#endif
