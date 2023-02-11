#ifndef __SHARED__
#define __SHARED__

#include <limits.h>
#include <time.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <sys/stat.h>
#include <float.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <math.h>

/* ========================================================================== */
/* MACRO DEFINITIONS ======================================================== */
/* ========================================================================== */

/* -> FUNCTION-LIKE MACROS ---------------------------------------------------*/

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))



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
 * except that DEBUG does not affect it (it is always active)
 * also prints the file and line info and exits the program
 * if the expression evaluates to false
 */
#define ASSERT(expr)                                                                             \
	if (!(expr))                                                                                 \
	{                                                                                            \
		fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); \
		exit(1);                                                                                 \
	}

/*
 * Macro:[ASSERTM]
 * shortcut to evaluate an expression, works the same way as ASSERT macro
 * but also prints a custom message
 */
#define ASSERTM(expr, msg)                                                                     \
	if (!(expr))                                                                               \
	{                                                                                          \
		fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, msg); \
		exit(1);                                                                               \
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

#define FREE2(expr, n)                   \
	if (expr)                            \
	{                                    \
		for (int i = 0; i < (int)n; i++) \
		{                                \
			free(expr[i]);               \
			expr[i] = NULL;              \
		}                                \
		free(expr);                      \
		expr = NULL;                     \
	}

/*
 * Macro:[FCLOSE]
 * shortcut to check if file is open and close it
 */
#define FCLOSE(expr)                                                                                 \
	if (expr)                                                                                        \
	{                                                                                                \
		if (fclose(expr) != 0)                                                                       \
		{                                                                                            \
			fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); \
			exit(1);                                                                                 \
		};                                                                                           \
		expr = NULL;                                                                                 \
	}

/*
 * Macro:[BGZCLOSE]
 * shortcut to check if bgzFile is open and close it
 */
#define BGZCLOSE(expr)                                                                               \
	if (expr)                                                                                        \
	{                                                                                                \
		if (bgzf_close(expr) != 0)                                                                   \
		{                                                                                            \
			fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); \
			exit(1);                                                                                 \
		};                                                                                           \
		expr = NULL;                                                                                 \
	}

/*
 * Macro:[GZCLOSE]
 * shortcut to check if gzFile is open and close it
 */
#define GZCLOSE(expr)                                                                                \
	if (expr)                                                                                        \
	{                                                                                                \
		if (gzclose(expr) != Z_OK)                                                                   \
		{                                                                                            \
			fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); \
			exit(1);                                                                                 \
		};                                                                                           \
		expr = Z_NULL;                                                                               \
	}

/* LIMIT DEFINING MACROS -----------------------------------------------------*/

// maximum number of tokens allowed in amova formula string
#define MAX_FORMULA_TOKENS 10

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

#define FREAD_BUF_SIZE 4096

#define FGETS_BUF_SIZE 4096

#define DELIMS "\t ,\n"

#define METADATA_DELIMS "\t ,\n"

#define MAX_N_AMOVA_LEVELS 5

const int pow10[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

const int nChoose2[301] = {0, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, 325, 351, 378, 406, 435, 465, 496, 528, 561, 595, 630, 666, 703, 741, 780, 820, 861, 903, 946, 990, 1035, 1081, 1128, 1176, 1225, 1275, 1326, 1378, 1431, 1485, 1540, 1596, 1653, 1711, 1770, 1830, 1891, 1953, 2016, 2080, 2145, 2211, 2278, 2346, 2415, 2485, 2556, 2628, 2701, 2775, 2850, 2926, 3003, 3081, 3160, 3240, 3321, 3403, 3486, 3570, 3655, 3741, 3828, 3916, 4005, 4095, 4186, 4278, 4371, 4465, 4560, 4656, 4753, 4851, 4950, 5050, 5151, 5253, 5356, 5460, 5565, 5671, 5778, 5886, 5995, 6105, 6216, 6328, 6441, 6555, 6670, 6786, 6903, 7021, 7140, 7260, 7381, 7503, 7626, 7750, 7875, 8001, 8128, 8256, 8385, 8515, 8646, 8778, 8911, 9045, 9180, 9316, 9453, 9591, 9730, 9870, 10011, 10153, 10296, 10440, 10585, 10731, 10878, 11026, 11175, 11325, 11476, 11628, 11781, 11935, 12090, 12246, 12403, 12561, 12720, 12880, 13041, 13203, 13366, 13530, 13695, 13861, 14028, 14196, 14365, 14535, 14706, 14878, 15051, 15225, 15400, 15576, 15753, 15931, 16110, 16290, 16471, 16653, 16836, 17020, 17205, 17391, 17578, 17766, 17955, 18145, 18336, 18528, 18721, 18915, 19110, 19306, 19503, 19701, 19900, 20100, 20301, 20503, 20706, 20910, 21115, 21321, 21528, 21736, 21945, 22155, 22366, 22578, 22791, 23005, 23220, 23436, 23653, 23871, 24090, 24310, 24531, 24753, 24976, 25200, 25425, 25651, 25878, 26106, 26335, 26565, 26796, 27028, 27261, 27495, 27730, 27966, 28203, 28441, 28680, 28920, 29161, 29403, 29646, 29890, 30135, 30381, 30628, 30876, 31125, 31375, 31626, 31878, 32131, 32385, 32640, 32896, 33153, 33411, 33670, 33930, 34191, 34453, 34716, 34980, 35245, 35511, 35778, 36046, 36315, 36585, 36856, 37128, 37401, 37675, 37950, 38226, 38503, 38781, 39060, 39340, 39621, 39903, 40186, 40470, 40755, 41041, 41328, 41616, 41905, 42195, 42486, 42778, 43071, 43365, 43660, 43956, 44253, 44551, 44850};

int find_n_given_nC2(int nC2_res);

/* ========================================================================== */
/* ENUMERATIONS ============================================================= */

// output file compression types
enum OUTFC
{
	NONE, // no compression [0]
	GZ,	  // gzip [1]
	BGZ,  // bgzip [2]
	BBGZ, // binary bgzip [3]
};

// input file types
enum INFT
{
	IN_VCF,
	IN_DM,
	IN_JGD,
	IN_JGPD,
};

/* ========================================================================== */

const double NEG_INF = -std::numeric_limits<double>::infinity();

int extractDigits(int num, int digits);

char *get_time();

void usage(FILE *fp);

#endif
