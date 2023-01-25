#ifndef __SHARED__
#define __SHARED__

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
#include <float.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <math.h>

// maximum number of tokens in amova formula string
#ifndef MAX_FORMULA_TOKENS
#define MAX_FORMULA_TOKENS 10
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/*
 * Macro:[MAXDIG_PER_HLEVEL]
 * Defines the maximum number of digits that can be used to
 * represent a single strata in a hierarchical level
 * e.g. MAX_DIGIT_PER_HLEVEL = 2
 *      strata0 = 00 (min)
 * 		...
 *      strata99 = 99 (max)
 *      nStrata = 100
 */
#define MAXDIG_PER_HLEVEL 2

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
#define FREE(expr)   \
	if (expr)        \
	{                \
		free(expr);  \
		expr = NULL; \
	}                \
	// fprintf(stderr, "\n\n*******\n[FREEING NULL MEMORY](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr);

/*
 * Macro:[FCLOSE]
 * shortcut to check if file is open and close it
 */
#define FCLOSE(expr)  \
	if (expr)         \
	{                 \
		fclose(expr); \
		expr = NULL;  \
	}                 \



#ifndef FREAD_BUF_SIZE
#define FREAD_BUF_SIZE 4096
#endif

#ifndef FGETS_BUF_SIZE
#define FGETS_BUF_SIZE 1024
#endif

#define DELIMS "\t ,\n"

#define METADATA_DELIMS "\t ,\n"

#define MAX_N_AMOVA_LEVELS 5

const int pow10[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

// input file type
enum
{
	IN_VCF,
	IN_DM,
	IN_SFS
};

const double NEG_INF = -std::numeric_limits<double>::infinity();

using size_t = decltype(sizeof(int));


int extractDigits(int num, int digits);


/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_vcf_fn		pointer to input file name
 * @field *in_mtd_fn	pointer to input Metadata file name
 * @field *out_fn		pointer to output file prefix [angsdput]
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
 * @field seed			random seed for bootstrapping
 * @field nBootstraps	number of bootstraps to be performed for AMOVA significance test
 *
 */
typedef struct
{

	int verbose;

	char *in_vcf_fn;
	char *in_sfs_fn;
	char *in_dm_fn;
	char *in_mtd_fn;
	char *out_fn;

	char *command;

	char *formula;
	int *keyCols;

	int blockSize;
	int doAMOVA;
	int doEM;

	int printMatrix;
	int isSim;
	int isTest;
	int doDist;
	int do_square_distance;
	int minInd;

	int hasColNames;

	double tole;
	int doTest;

	int mThreads;
	int mEmIter;

	int seed;
	int nBootstraps;

	int gl2gt;

} argStruct;

argStruct *argStruct_init();
argStruct *argStruct_get(int argc, char **argv);
void argStruct_destroy(argStruct *arg);
void argStruct_print(FILE* fp, argStruct *arg);

/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites				number of sites
 * @field nInd					number of individuals
 * @field pos					position
 *
 * @field LUT_indPairIdx		lookup table for mapping two individuals
 * 								to their pair index
 * @field nIndCmb				number of unique pairwise individual combinations
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

	// int *pos;

	int **LUT_indPairIdx;
	int nIndCmb;

	// char *major;
	// char *minor;
	// char *ref;
	// char *anc;
	// char *der;

	// input file type from enum
	int in_ft;

	char *DATETIME;

	void print()
	{
		// print lookup table
		for (int i1 = 0; i1 < nInd - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < nInd; i2++)
			{
				fprintf(stderr, "\n%i %i %i\n", LUT_indPairIdx[i1][i2], i1, i2);
			}
		}
	}

	// validate that parameters make sense
	void validate()
	{
		ASSERT(nIndCmb > 0);
		ASSERT(nInd > 0);
		ASSERT(nSites > 0);
		ASSERT(totSites > 0);
	}

} paramStruct;

paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

char *get_time();

void usage(FILE *fp);

/**
 * @namespace DATA
 * @abstract data structures
 *
 * @typedef formulaStruct - structure for storing AMOVA formula
 *
 */
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
				fprintf(fp, "\n\t%i\t%s\n", i, formulaTokens[i]);
			}
		}

	} formulaStruct;
	formulaStruct *formulaStruct_get(const char *formula);


	/// @brief blobStruct - structure for storing blocks of contig data
	typedef struct blobStruct
	{

		// Number of contigs
		size_t nContigs;

		// array of contig indices
		int *contigIdx;

		// 2d array of contig names
		// [Contig][Name]
		char **contigNames;

		// Array of number of blocks per contig
		// contigNBlocs[Contig] = number of blocks in the contig
		int *contigNBlocks;

		// Array of contig lengths
		// contigLengths[Contig] = length of the contig
		int *contigLengths;

		// 2D array of contig block starts
		// [Contig][BlockStart]
		int **contigBlockStarts;

		// [Contig][Block] int* to block start in vcfData->lngl
		// 2D array of pointers to actual contig block starts
		double ***contigBlockStartPtrs;


	} blobStruct;

	blobStruct *blobStruct_init(const int nContigs, const int blockSize, bcf_hdr_t *hdr);
	void blobStruct_destroy(blobStruct *c);

	typedef struct pairStruct
	{

		int i1;
		int i2;
		int idx;

		// number of shared sites
		size_t snSites = 0;

		// contains index of the shared sites of the pair
		int *sharedSites = NULL;
		size_t _sharedSites = 1024;


		double d;
		int n_em_iter;

		// TODO should not always create below
		//  double SFS[3][3]={{NEG_INF,NEG_INF,NEG_INF},
		//  {NEG_INF,NEG_INF,NEG_INF},
		//  {NEG_INF,NEG_INF,NEG_INF}};

		double *SFS = NULL;

		pairStruct()
		{
			d = 0.0;
			n_em_iter = 0;
			i1=-1;
			i2=-1;
			idx=-1;
			sharedSites = (int *)malloc(_sharedSites * sizeof(int));
			for (size_t i = 0; i < _sharedSites; i++)
			{
				sharedSites[i] = -1;
			}

			SFS = (double *)malloc(9 * sizeof(double));
			for (int i = 0; i < 9; i++)
			{
				SFS[i] = NEG_INF;
			}
		}
		~pairStruct()
		{
			FREE(sharedSites);
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

		void sharedSites_add(size_t site_i)
		{
			sharedSites[snSites] = site_i;
			snSites++;
			if (snSites == _sharedSites)
			{
				sharedSites_expand();
			}
		}

		void sharedSites_expand()
		{
			_sharedSites *= 2;
			sharedSites = (int *)realloc(sharedSites, _sharedSites * sizeof(int));
			for (size_t i = snSites; i < _sharedSites; i++)
			{
				sharedSites[i] = -1;
			}
		}

	} pairStruct;



	/// @brief sampleStruct - contains sample names and other sample related information
	/// used in comparing samples in the VCF file with the samples in the metadata file
	typedef struct sampleStruct
	{

		char **sampleNames = NULL; // sample names in bcf sample order

		int nSamples;

		sampleStruct()
		{
			nSamples=-1;
		}

		~sampleStruct()
		{
			for (int i = 0; i < nSamples; i++)
			{
				free(sampleNames[i]);
				sampleNames[i] = NULL;
			}
			free(sampleNames);
			sampleNames = NULL;
		}

		void init(int nSamples_)
		{
			nSamples = nSamples_;

			sampleNames = (char **)malloc(nSamples * sizeof(char *));

			for (size_t i = 0; i < (size_t)nSamples; i++)
			{
				sampleNames[i] = NULL;
			}
		}

		void addSampleName(int i, char *str)
		{
			size_t size = strlen(str) + 1;
			sampleNames[i] = (char *)malloc(size);
			strncpy(sampleNames[i], str, size);
		}

	} sampleStruct;

	//TODO calculate associations once and store in a LUT
	/// @brief hierStruct store the hierarchical structure of the metadata
	typedef struct hierStruct
	{

		// number of unique stratas at this level
		// e.g. Population = {POP1, POP2, POP3}
		// 		nStrata = 3
		int nStrata = 0;

		// strata names
		char **strataNames = NULL;
		size_t _strataNames = 1;

		char *strataArr = NULL;

		int *nIndPerStrata = NULL;
		size_t _nIndPerStrata = 1024;

		// hierStruct(char* str, int totNumLevels)
		hierStruct(char *str)
		{
			nStrata = 1;
			size_t size = strlen(str) + 1;
			strataNames = (char **)malloc(_strataNames * sizeof(char *));
			strataNames[0] = (char *)malloc(size);
			strncpy(strataNames[0], str, size);

			nIndPerStrata = new int[_nIndPerStrata]; // TODO
			for (size_t i = 0; i < _nIndPerStrata; i++)
			{
				nIndPerStrata[i] = 0;
			}

			// nMemberStrata = new int*[totNumLevels];
			// for (size_t i = 0; i < totNumLevels; i++)
			// {
			// 	nMemberStrata[i] = NULL;
			// }
		}

		~hierStruct()
		{
			for (size_t i = 0; i < _strataNames; i++)
			{
				free(strataNames[i]);
				strataNames[i] = NULL;
			}
			free(strataNames);
			strataNames = NULL;

			delete[] nIndPerStrata;
		}

		void addStrata(char *str)
		{

			_strataNames = nStrata + 1;

			strataNames = (char **)realloc(strataNames, _strataNames * sizeof(char *));
			ASSERT(strataNames != NULL);

			strataNames[nStrata] = (char *)malloc(strlen(str) + 1);
			ASSERT(strataNames[nStrata] != NULL);

			strncpy(strataNames[nStrata], str, strlen(str) + 1);
			nIndPerStrata[nStrata] = 1;

			nStrata++;
		}

		int getStrataIndex(char *str)
		{
			int idx = 0;
			// check the current records
			while (idx < nStrata)
			{
				if (strcmp(str, strataNames[idx]) == 0)
				{
					nIndPerStrata[idx]++;
					return idx;
				}
				idx++;
			}

			// if not found
			addStrata(str);
			return idx;
		}

	} hierStruct;

	typedef struct metadataStruct
	{

		// number of hierarchical levels
		// e.g. Individual, Region, Population, Subpopulation
		// 		has 3 levels: {Individual, Region, Population}
		int nLevels = 0;

		char **levelNames = NULL;
		size_t _levelNames = 0;

		// associate individual with the lowest level strata (key strata)
		size_t *ind2stratakey = NULL;
		size_t _ind2stratakey = 1024; // TODO

		int nIndMetadata = 0;

		// size_t *strataKeys = NULL;

		// hierStruct[hierarchyLevel] is a hierStruct instance
		// exclude the lowest level (Individual)
		// Individual ~ Region / Population
		// hierArr = {Region, Population}
		hierStruct **hierArr = NULL;

		metadataStruct(int nLevels_, char **levelNames_)
		{
			nLevels = nLevels_;
			_levelNames = nLevels + 1;

			hierArr = new hierStruct *[nLevels];
			for (int i = 0; i < nLevels; i++)
			{
				hierArr[i] = NULL;
			}

			levelNames = new char *[_levelNames];
			for (size_t i = 0; i < _levelNames; i++)
			{
				levelNames[i] = new char[strlen(levelNames_[i]) + 1];
				strncpy(levelNames[i], levelNames_[i], strlen(levelNames_[i]) + 1);
			}

			ind2stratakey = new size_t[_ind2stratakey]; // TODO nInd

			// maximum number of strata keys allowed
			// MAXDIG_PER_HLEVEL = 2 (default) -> size_t[99]
			// strataKeys = new size_t[pow(10,MAXDIG_PER_HLEVEL)-1];
		};

		~metadataStruct()
		{
			for (int i = 0; i < nLevels; i++)
			{
				delete hierArr[i];
				delete[] levelNames[i];
			}
			delete[] levelNames[nLevels]; // levelNames if of size nLevels+1
			delete[] levelNames;
			delete[] hierArr;
			delete[] ind2stratakey;
		};

		void addHierStruct(size_t lvl_i, char *tok)
		{
			hierArr[lvl_i] = new hierStruct(tok);
		}

		// check if two individuals are in the same strata at a given level
		int pairInSameStrataAtLevel(int lvl, int i1, int i2)
		{
			if (getDigits(ind2stratakey[i1], lvl) == getDigits(ind2stratakey[i2], lvl))
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}

		void printLevelNames()
		{
			for (int i = 0; i < nLevels + 1; i++)
			{
				fprintf(stderr, "\n\n[INFO]\t-> levelNames[%d] = %s\n", i, levelNames[i]);
			}
		}

		/// @brief calculateKeyAtLevel - calculate the strata key value at a given strata index in a given level
		/// @param lvl			- hierarchical level
		/// @param strata_idx 	- strata index
		/// @return integer value of the strata key
		/// @example lvl=1, strata_idx=2 -> 1[02]00
		int calculateKeyAtLevel(int lvl, int strata_idx)
		{
			return (int)(pow(10, MAXDIG_PER_HLEVEL * (lvl +1)) + strata_idx);
			// pow10[]
		}

		int pairInStrataAtLevel(int i1, int i2, int lvl, int strata_idx)
		{
			int digit_idx = nLevels - lvl -1;
			if( (strata_idx==getDigits(ind2stratakey[i1], digit_idx)) && (strata_idx==getDigits(ind2stratakey[i2], digit_idx)) )
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}

		void print(FILE *fp)
		{
			fprintf(fp, "\nnLevels: %d\n", nLevels);
			for (int i = 0; i < nLevels; i++)
			{
				// fprintf(fp, "Level %d: contains %d unique group identifiers\n", i, hierArr[i]->nStrata);
				fprintf(fp, "Level (index:%d,id:%s) contains %d unique group identifiers: {", i, levelNames[i], hierArr[i]->nStrata);
				for (int j = 0; j < hierArr[i]->nStrata; j++)
				{
					fprintf(fp, "%s", hierArr[i]->strataNames[j]);
					if (j < hierArr[i]->nStrata - 1)
					{
						fprintf(fp, ",");
					}
					else
					{
						fprintf(fp, "}\n");
					}
				}
				for (int j = 0; j < hierArr[i]->nStrata; j++)
				{
					fprintf(fp, "\tGroup (index:%d,id:%s) has %d members\n", j, hierArr[i]->strataNames[j], hierArr[i]->nIndPerStrata[j]);
					for (int k = 0; k < hierArr[i]->nIndPerStrata[j]; k++)
					{
						fprintf(fp, "\t\tMember %d belongs to %s\n", k, hierArr[i]->strataNames[j]);
					}
				}
			}
		}

		void print_ind2stratakey(FILE *fp)
		{
			fprintf(fp, "\n\n------------------\n");
			fprintf(stderr, "\n\n nIndMetadata: %d\n", nIndMetadata);
			for (int i = 0; i < nIndMetadata; i++)
			{
				fprintf(fp, "ind=%d:stratakey=%ld\n", i, ind2stratakey[i]);
			}
			fprintf(fp, "\n------------------\n\n");
		};

		void print_pair2assoc(FILE *fp){

			for(int lvl=0; lvl < nLevels; lvl++){
				for (int sti=0; sti < hierArr[lvl]->nStrata; sti++)
				{

					for (int i1 = 0; i1 < nIndMetadata - 1; i1++)
					{
						for (int i2 = i1 + 1; i2 < nIndMetadata; i2++)
						{
							// include only pairs that are in the same strata at this hierarchical level
							if(pairInStrataAtLevel(i1, i2, lvl, sti) == 1){
							fprintf(stderr, "\n->Pair %i %i ", i1, i2);
								fprintf(stderr, "in strata %i at level %d", sti,lvl);
							}
						}
					}
				}
			}

		}

		/// @brief get the number at a given digit of a number
		/// @param number	number to extract digit from
		/// @param digit	digit to extract (right to left, 0-based)
		// digits [digit_2][digit_1][digit_0]
		//   	  100 -> digit_2 = 1, digit_1 = 0, digit_0 = 0
		/// @return extracted digit
		int getDigit(int number, int digit)
		{
			if (number < (int)pow(10, digit))
			{
				return 0;
			}
			else
			{
				return (number / (int)pow(10.0, digit)) % 10;
			}
		}

		/// @brief getDigits - extract n_digits consecutive digits from number starting at idx (right to left)
		/// @param number	number to extract digits from
		/// @param idx		index of the first digit to extract
		/// @param n_digits	number of consecutive digits to extract (starting at idx, right to left)
		/// @return extracted digits
		/// @example number=99123456, idx=1, n_digits=2 -> [99][12][34][56] -> 34
		int getDigits(int number, int idx, int n_digits)
		{


			int res = number % pow10[(n_digits* idx)+n_digits];

			if(idx>0)
			{
				res = res / pow10[(n_digits* (idx-1))+n_digits];
			}

			return res;
		}

		/// @brief getDigits - overload default value for n_digits
		/// @param number	number to extract digits from
		/// @param idx		index of the first digit to extract
		/// @return
		int getDigits(int number, int idx)
		{

			int res = number % pow10[(MAXDIG_PER_HLEVEL* idx)+2];

			if(idx>0)
			{
				res = res / pow10[(MAXDIG_PER_HLEVEL* (idx-1))+2];
			}
			return res;
		}


		size_t initKey()
		{
			size_t key = (size_t)(1 * pow(10, MAXDIG_PER_HLEVEL * (nLevels)));
			return key;
		}

		size_t setKeyDigitAtLevel(size_t key, int lvl_i, int lvl_strata_i)
		{

			if (MAXDIG_PER_HLEVEL != 2)
			{
				fprintf(stderr, "\n[ERROR] MAXDIG_PER_HLEVEL currently only supports 2\n");
				exit(1);
			}

			int dig_lvl= nLevels - lvl_i -1;

			size_t ret = pow10[dig_lvl * MAXDIG_PER_HLEVEL] * lvl_strata_i;
			ret = ret + key;

			return ret;
		}

		// 123456789 lvl=2 MAXDIG_PER_HLEVEL=2
		// get lvl 2 and lower hier levels
		// == 456789
		size_t extractSubkey(size_t key, int lvl)
		{

			size_t ret = key % pow10[(MAXDIG_PER_HLEVEL * (lvl) + 2)];
			// if(ret == 0){
				// if 910023 % 10000 = 0, return 0 + 10000
				ret = ret + pow10[(MAXDIG_PER_HLEVEL * (lvl) + 2)];
			// }
			return ret;
		}

		// TODO
		/// @brief keyIsFromSameStrataGivenStrataLevel - check if key is from the same strata at level
		/// @param key
		/// @param strata_i
		/// @param lvl
		/// @return  1 if key is from the same strata at level, 0 otherwise
		/// @example 12345 lvl=1 strata_i=23
		///                calculateKeyAtLevel(strata_i=23,lvl=1) = 2300
		/// 			   extractSubkey(key=12345,lvl=1) = 2345
		/// 			   2299 < 2345 < 2400
		///
		int keyIsFromSameStrataGivenStrataLevel(size_t key, int strata_i, int lvl)
		{
			int k_i = extractSubkey(key, lvl);
			int ref_k = calculateKeyAtLevel(lvl, strata_i);
			// assume MAXDIG_PER_HLEVEL = 2
			if (ref_k - 1 < k_i && k_i < ref_k + 100)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}

		int keyIsFromSameStrataAsKey(size_t key, size_t ref_key, int lvl)
		{
			int k_i = extractSubkey(key, lvl);
			int ref_k = extractSubkey(ref_key, lvl);
			// assume MAXDIG_PER_HLEVEL = 2
			if (ref_k - 1 < k_i && k_i < ref_k + 100)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}

		int nMemberStrataAtStrataAtLevel(int strata_i, int lvl)
		{
			//digits count from right side
			// hlvl0     hlvl1     hlvl2
			// digitlvl2 digitlvl1 digitlvl0
			int diglvl= nLevels-1-lvl;
			
			ASSERT(diglvl>0);
			size_t key_i = 0;

			int found = 0;
			int sum = 0;
			int tmp[pow10[MAXDIG_PER_HLEVEL]];

			// loop through existing keys (one key per individual)
			for (int i = 0; i < nIndMetadata; i++)
			{
				key_i = ind2stratakey[i];

				if(getDigits(key_i,diglvl) == strata_i){
					int new_subkey = extractSubkey(key_i, diglvl - 1);

					found = 0;
					for (int j = 0; j < sum; j++)
					{
						if (tmp[j] == new_subkey)
						{
							found = 1;
							break;
						}
					}
					if (found == 1)
					{
						continue;
					}
					else
					{
						tmp[sum] = new_subkey;
						sum++;
					}
				}
			}
			return sum;
		}


		//
		// for each strata at this level;
		// 		number of unique member stratas from lowel levels
		//
		// Metadata:
		// Individual, Region, Population, Subpopulation
		// ind1, reg1, pop1, subpop1
		// ind2, reg1, pop1, subpop1
		// ind3, reg1, pop1, subpop2
		// ind4, reg1, pop2, subpop3
		// ind5, reg1, pop2, subpop3
		// ind6, reg2, pop3, subpop4
		// ind7, reg2, pop3, subpop4
		//
		// hierStruct[0] = Region
		// hierStruct[0]->nMemberStrata = {memberPopulation, memberSubpopulation}
		// hierStruct[0]->nMemberStrata[0] = memberPopulation = {2,1}
		//    2 member stratas {pop1,pop2} for reg1
		//    1 member strata {pop3} for reg2
		//
		// hierStruct[0]->nMemberStrata[1] = memberSubpopulation = {3,1}
		//   3 member stratas {subpop1,subpop2,subpop3} for reg1
		//   1 member strata {subpop4} for reg2
		//
		// nMemberStrata[Number of hlevels lower than this hlevel]
		int sumUniqStrataAtEachLevel(int lvl)
		{
			if (lvl == nLevels - 1)
			{
				return hierArr[0]->nStrata;
			}else if (nLevels == 1)
			{
				ASSERT(0 == 1);
			}else if (nLevels > 1 && lvl < nLevels - 1)
			{

				int lvl_i = lvl + 1;
				int sum = 0;
				for (int i = 0; i < hierArr[lvl_i]->nStrata; i++)
				{

					sum += nMemberStrataAtStrataAtLevel(i, lvl_i-1);
				}
				return sum;
			}
			ASSERT(0 == 1);
		}

	} metadataStruct;

	metadataStruct *metadataStruct_get(FILE *in_mtd_fp, sampleStruct *SAMPLES, formulaStruct *FORMULA, int has_colnames, paramStruct *pars);

	/**
	 * @brief distanceMatrixStruct stores the distance matrix
	 *
	 * @param nInd 		number of individuals
	 * @param nIndCmb	number of individual combinations
	 * @param isSquared 1 if the distance matrix is squared, 0 otherwise
	 *
	 */
	typedef struct distanceMatrixStruct
	{
		double *M = NULL;


		int nInd = 0;
		int nIndCmb = 0;
		int isSquared = -1;

		distanceMatrixStruct(int nInd_, int nIndCmb_, int isSquared_)
		{
			nIndCmb = nIndCmb_;
			nInd = nInd_;
			M = new double[nIndCmb];
			for (int i = 0; i < nIndCmb; i++)
			{
				M[i] = 0;
			}
			isSquared = isSquared_;
		};
		~distanceMatrixStruct()
		{
			delete[] M;
		}
		void print(FILE *fp)
		{
			for (int px = 0; px < nIndCmb; px++)
			{
				fprintf(fp, "%.*f", (int)DBL_MAXDIG10, M[px]);
				if (px != nIndCmb - 1)
				{
					fprintf(fp, ",");
				}
				else
				{
					fprintf(fp, "\n");
				}
			}
		}

		void print(FILE *fp, const char *TYPE)
		{
			fprintf(fp, "%s,", TYPE);
			for (int px = 0; px < nIndCmb; px++)
			{
				fprintf(fp, "%.*f", (int)DBL_MAXDIG10, M[px]);
				if (px != nIndCmb - 1)
				{
					fprintf(fp, ",");
				}
				else
				{
					fprintf(fp, "\n");
				}
			}
		}
	} distanceMatrixStruct;

	// read distance matrix from distance matrix file
	distanceMatrixStruct *distanceMatrixStruct_read(FILE *in_dm_fp, paramStruct *pars, argStruct *args);
};

namespace IO
{

	char *setFileName(const char *a, const char *b);

	FILE *getFile(const char *fname, const char *mode);
	FILE *openFile(const char *a, const char *b);
	FILE *openFile(char *c);

	namespace readFile
	{
		int SFS(FILE *in_sfs_fp, const char *delims, DATA::sampleStruct *SAMPLES);

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
		FILE *fp = NULL;
		// int has_header = 0; //TODO

		outputStruct(const char *fn_, const char *suffix)
		{
			fn = setFileName(fn_, suffix);
			fp = openFile(fn);
		}
		~outputStruct()
		{
			fclose(fp);
			free(fn);
			fn = NULL;
		}
		void flush()
		{
			fflush(fp);
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

			if (args->doEM == 1)
			{
				out_sfs_fs = new outputStruct(args->out_fn, ".sfs.csv");
			}
			if (args->doAMOVA > 0)
			{
				out_amova_fs = new outputStruct(args->out_fn, ".amova.csv");
			}
		}

		~outFilesStruct()
		{
			delete out_emtest_fs;
			delete out_dm_fs;
			delete out_sfs_fs;
			delete out_amova_fs;
		}

		void flushAll()
		{
			if (out_emtest_fs != NULL)
			{
				out_emtest_fs->flush();
			}
			if (out_dm_fs != NULL)
			{
				out_dm_fs->flush();
			}
			if (out_sfs_fs != NULL)
			{
				out_sfs_fs->flush();
			}
			if (out_amova_fs != NULL)
			{
				out_amova_fs->flush();
			}
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

		void M_PWD(const char *TYPE, IO::outputStruct *out_dm_fs, int nIndCmb, double *M_PWD);
	}
}

typedef struct threadStruct
{

	DATA::pairStruct *pair;
	double **lngls;

	FILE *out_sfs_fp;

	// TODO use them globally?
	double tole;
	int mEmIter;

	size_t nSites;

	threadStruct(DATA::pairStruct *tPair, double **lngl, size_t nSites_t, IO::outputStruct *out_sfs_fs, argStruct *args)
	{
		pair = tPair;
		lngls = lngl;
		nSites = nSites_t;
		out_sfs_fp = out_sfs_fs->fp;
		tole = args->tole;
		mEmIter = args->mEmIter;
	}

} threadStruct;
void print_SFS_GT(const char *TYPE, IO::outputStruct *out_sfs_fs, paramStruct *pars, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);


void check_consistency_args_pars(argStruct *args, paramStruct *pars);

#endif
