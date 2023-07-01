#ifndef __SHARED__
#define __SHARED__

#include <ctype.h>
#include <float.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include <limits>

#include "argStruct.h"
#include "dev.h"
#include "lookup.h"

/* ========================================================================== */
/* MACRO DEFINITIONS ======================================================== */
/* ========================================================================== */

/* -> FUNCTION-LIKE MACROS ---------------------------------------------------*/

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/*
 * Macro:[DBL_MAX_DIG_TOPRINT]
 * 	maximum number of digits needed to print a double
 *
 * 	longest number == smalles negative number
 * 		-pow(2, DBL_MIN_EXP - DBL_MANT_DIG)
 * 	-pow(2,-N) needs 3+N digits
 * 		to represent (sign, decimal point, N digits)
 * 		'-0.<N digits>'
 *
 * @requires <float.h>
 */
#define DBL_MAX_DIG_TOPRINT 3 + DBL_MANT_DIG - DBL_MIN_EXP
// TODO deprecated

/*
 * Macro:[AT]
 * inject the file and line info as string
 */
#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)

/*
 * Macro:[ERROR]
 * print a custom error message and exit the program
 */
#define ERROR(...)                                                                           \
    fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
    fprintf(stderr, __VA_ARGS__);                                                            \
    fprintf(stderr, "\n*******\n");                                                          \
    exit(1);

/*
 * Macro:[NEVER]
 * indicates that a point in the code should never be reached
 */
#define NEVER \
    ERROR("Control should never reach this point; please report this to the developers.")

/*
 * Macro:[ASSERT]
 * evaluate an expression, works the same way as the C-macro assert
 * except that DEBUG does not affect it (it is always active)
 * also prints the file and line info and exits the program
 * if the expression evaluates to false
 */
#define ASSERT(expr)                                                                                              \
    if (!((expr))) {                                                                                              \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", __FILE__, __FUNCTION__, __LINE__, #expr); \
        exit(1);                                                                                                  \
    }

/*
 * Macro:[WARNING]
 * print a custom warning message
 */
#define WARNING(...)                                                                \
    fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, __FUNCTION__, __LINE__); \
    fprintf(stderr, __VA_ARGS__);                                                   \
    fprintf(stderr, "\n");                                                          \
    IO::append_argsFile("\n[WARNING](%s/%s:%d): %s\n", __FILE__, __FUNCTION__, __LINE__, __VA_ARGS__);

/*
 * Macro:[LOG]
 */
#define LOG(...)                                    \
    fprintf(stderr, "\n[LOG](%s)\t", __FUNCTION__); \
    fprintf(stderr, __VA_ARGS__);

#define BEGIN_LOGSECTION                                             \
    fprintf(stderr, "\n________________________________________\n"); \
    fprintf(stderr, "[LOG]\t-> Program is now at <%s/%s>\n", __FILE__, __FUNCTION__);

#define END_LOGSECTION \
    fprintf(stderr, "\n________________________________________\n");

#define ARGLOG(...)                                 \
    fprintf(stderr, "\n[LOG](%s)\t", __FUNCTION__); \
    fprintf(stderr, __VA_ARGS__);                   \
    IO::append_argsFile("\n[LOG](%s): %s\n", __FUNCTION__, __VA_ARGS__);

/*
 * Macro:[ASSERTM]
 * shortcut to evaluate an expression, works the same way as ASSERT macro
 * but also prints a custom message
 */
#define ASSERTM(expr, msg)                                                                               \
    if (!(expr)) {                                                                                       \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n", __FILE__, __FUNCTION__, __LINE__, #expr); \
        fprintf(stderr, "\n%s!!\n", #msg);                                                               \
        exit(1);                                                                                         \
    }

/*
 * Macro:[CREALLOC] - check if realloc was successful
 *
 * if so, set original pointer to new memory
 * otherwise print error message and exit
 *
 * e.g. int **rc_p = realloc(p, sizeof(int) * 10);
 *     CREALLOC(p);
 *
 */
#define CREALLOC(p)                                                                              \
    if (NULL == rc_##p) {                                                                        \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
        fprintf(stderr, "Failed to reallocate memory\n");                                        \
        fprintf(stderr, "\n*******\n");                                                          \
        exit(1);                                                                                 \
    } else {                                                                                     \
        p = rc_##p;                                                                              \
    }

/*
 * Macro:[CCREALLOC] - check if realloc was successful
 *
 * if so, set original pointer to new memory
 * otherwise print error message and exit
 *
 * e.g. int **r_p = realloc(p, sizeof(int) * 10);
 *     CREALLOC(r_p, p);
 *
 *
 */
#define CCREALLOC(rp, p)                                                                         \
    if (NULL == rp) {                                                                            \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
        fprintf(stderr, "Failed to reallocate memory\n");                                        \
        fprintf(stderr, "\n*******\n");                                                          \
        exit(1);                                                                                 \
    } else {                                                                                     \
        p = rp;                                                                                  \
    }

/*
 * Macro:[FREE]
 * force free
 *
 * shortcut to free memory and set pointer to NULL
 * and catch double free
 * throw error if pointer is NULL
 */
#define FREE(x)                                 \
    if (NULL != x) {                            \
        free(x);                                \
        x = NULL;                               \
    } else {                                    \
        ERROR("Trying to free a NULL pointer"); \
    }

/*
 * Macro:[FFREE]
 * if x exists; free memory and set pointer to NULL
 * do NOT throw error if pointer is NULL
 */
#define FFREE(x)     \
    if (NULL != x) { \
        free(x);     \
        x = NULL;    \
    }

/*
 * Macro:[FREE2D]
 * free 2D array
 * throw error if pointer is NULL
 */
#define FREE2D(x, size)                                 \
    if (NULL != x) {                                    \
        for (size_t i = 0; i < (size_t)size; i++) {     \
            if (NULL != x[i]) {                         \
                free(x[i]);                             \
                x[i] = NULL;                            \
            } else {                                    \
                ERROR("Trying to free a NULL pointer"); \
            }                                           \
        }                                               \
        free(x);                                        \
        x = NULL;                                       \
    } else {                                            \
        ERROR("Trying to free a NULL pointer");         \
    }

/*
 * Macro:[IFFREE2D]
 * if x exists; free 2D array
 * do NOT throw error if pointer is NULL
 */
#define IFFREE2D(x, size)                           \
    if (NULL != x) {                                \
        for (size_t i = 0; i < (size_t)size; i++) { \
            if (NULL != x[i]) {                     \
                free(x[i]);                         \
                x[i] = NULL;                        \
            }                                       \
        }                                           \
        free(x);                                    \
        x = NULL;                                   \
    }

/*
 * Macro:[DEL]
 * delete memory and set pointer to NULL
 * throw error if pointer is NULL
 *
 * for variables declared as:
 * type *x = new type;
 */
#define DEL(x)                                    \
    if (NULL != x) {                              \
        delete x;                                 \
        x = NULL;                                 \
    } else {                                      \
        ERROR("Trying to delete a NULL pointer"); \
    }

/*
 * Macro:[FDEL]
 * if x exists; delete memory and set pointer to NULL
 * do NOT throw error if pointer is NULL
 *
 * for variables declared as:
 * type *x = new type;
 */
#define FDEL(x)      \
    if (NULL != x) { \
        delete x;    \
        x = NULL;    \
    }

/*
 * Macro:[DEL1D]
 * delete memory and set pointer to NULL
 * throw error if pointer is NULL
 *
 * for 1D arrays declared as:
 * type *x = new type[size];
 *
 */
#define DEL1D(x)                                  \
    if (NULL != x) {                              \
        delete[] x;                               \
        x = NULL;                                 \
    } else {                                      \
        ERROR("Trying to delete a NULL pointer"); \
    }

/*
 * Macro:[FDEL1D]
 * if x exists; delete memory and set pointer to NULL
 * do NOT throw error if pointer is NULL
 *
 * for 1D arrays declared as:
 * type *x = new type[size];
 *
 */
#define FDEL1D(x)    \
    if (NULL != x) { \
        delete[] x;  \
        x = NULL;    \
    }

/*
 * Macro:[DEL2D]
 * shortcut to delete memory and set pointer to NULL
 * for 2D arrays declared as:
 * type **x = new type*[size];
 * for (size_t i = 0; i < size; i++) {
 *   x[i] = new type[size];
 * }
 */
#define DEL2D(x, size)                              \
    if (NULL != x) {                                \
        for (size_t i = 0; i < (size_t)size; i++) { \
            if (NULL != x[i]) {                     \
                delete x[i];                        \
                x[i] = NULL;                        \
            }                                       \
        }                                           \
        delete[] x;                                 \
        x = NULL;                                   \
    } else {                                        \
        ERROR("Trying to delete a NULL pointer");   \
    }

/*
 * Macro:[FDEL2D]
 * shortcut to delete memory and set pointer to NULL
 * for 2D arrays declared as:
 * type **x = new type*[size];
 * for (size_t i = 0; i < size; i++) {
 *   x[i] = new type[size];
 * }
 *
 * do NOT throw error if pointer is NULL
 */
#define FDEL2D(x, size)                             \
    if (NULL != x) {                                \
        for (size_t i = 0; i < (size_t)size; i++) { \
            if (NULL != x[i]) {                     \
                delete x[i];                        \
                x[i] = NULL;                        \
            }                                       \
        }                                           \
        delete[] x;                                 \
        x = NULL;                                   \
    }

/*
 * Macro:[FCLOSE]
 * shortcut to check if file is open and close it
 */
#define FCLOSE(expr)                                                                                 \
    if (expr) {                                                                                      \
        if (fclose(expr) != 0) {                                                                     \
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
    if (expr) {                                                                                      \
        if (bgzf_close(expr) != 0) {                                                                 \
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
    if (expr) {                                                                                      \
        if (gzclose(expr) != Z_OK) {                                                                 \
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #expr); \
            exit(1);                                                                                 \
        };                                                                                           \
        expr = Z_NULL;                                                                               \
    }

/* -> BIT MANIPULATION MACROS ------------------------------------------------*/

/* Macro:[BITSET]
 * set a specific bit in x
 *
 * @param x		the variable to set the bit in
 * @param bit	the bit to set
 * @return		x is modified in place
 */
// #define BITSET(x, bit) ((x) |= (1ULL << (bit)))
#define BITSET(x, bit)                                      \
    if ((bit) < 0 || (bit) >= 64) {                         \
        fprintf(stderr, "Invalid bit number: %d\n", (bit)); \
        exit(1);                                            \
    }                                                       \
    (x) |= (1ULL << (bit));

/* Macro:[BITTOGGLE]
 * toggle a specific bit in x
 *
 * @param x		the variable to toggle the bit in
 * @param bit	the bit to toggle
 * @return		x is modified in place
 */
#define BITTOGGLE(x, bit) ((x) ^= (1ULL << (bit)))

/* Macro:[BITCLEAR]
 * clear a specific bit in x
 *
 * @param x		the variable to clear the bit in
 * @param bit	the bit to clear
 * @return		x is modified in place
 */
#define BITCLEAR(x, bit) ((x) &= ~(1ULL << (bit)))

/* Macro:[BITCHECK]
 * check if a specific bit is set in x
 *
 * @param x		the variable to check the bit in
 * @param bit	the bit to check
 * @return		1 if the bit is set, 0 otherwise
 */
#define BITCHECK(x, bit) !!((x) & (1ULL << (bit)))

/* Macro:[CHAR_BITCHECK_ANY]
 * check if any bit is set in a char
 *
 * @param x		the char to check
 * @return		1 if any bit is set, 0 otherwise
 */
#define CHAR_BITCHECK_ANY(x) !!((x)&0xFF)

/* Macro:[WHICH_BIT_SET]
 *
 * @param x		the variable to check
 * @return		the index of the first bit set in x
 * 				-1 if no bit is set
 */
#define WHICH_BIT_SET(x) (x == 0 ? -1 : (int)log2(x))

/* Macro:[WHICH_BIT_SET1]
 *
 * @param x		the variable to check
 * @return		the 1-based index of the first bit set in x
 * 				-1 if no bit is set
 */
#define WHICH_BIT_SET1(x) (x == 0 ? -1 : (int)log2(x) + 1)

/* Macro:[BITCHECK_ATLEAST]
 * check if at any bit that is at least as significant as the specified bit is set
 *
 * @param x		the variable to check
 * @param bit	the bit to check
 * @return		1 if any bit that is at least as significant as the specified bit is set
 * 				0 otherwise
 */
#define BITCHECK_ATLEAST(x, bit) !!((x >> (bit)) & 0xFF)

/* LIMIT DEFINING MACROS -----------------------------------------------------*/

// maximum number of tokens allowed in amova formula string
#define MAX_N_FORMULA_TOKENS 10

#define MAX_N_HIER_LEVELS 10
// TODO
#define MAX_N_GROUPS_PER_LEVEL 10

#define MAX_N_BITS 64

// maximum level of verbosity allowed
#define MAX_VERBOSE_LEVEL 7

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
#define DBL_MAXDIG10 ((int)(2 + (DBL_MANT_DIG * 30103UL) / 100000UL))

#define METADATA_DELIMS "\t ,\n"

// dragon
#define BUF_NSITES 4096
#define BUF_NCONTIGS 256
#define BUF_NTOTSITES 4096

#define FREAD_BUF_SIZE 4096
#define FGETS_BUF_SIZE 25600

/* CONSTANTS -----------------------------------------------------------------*/

// 1/9
#define FRAC_1_9 0.1111111111111111

/* ========================================================================== */
/* ENUMERATIONS ============================================================= */

// output file compression types
enum OUTFC {
    NONE,  // no compression [0]
    GZ,    // gzip [1]
    BBGZ,  // bgzip (binary) [2]
};

// input file types for main data input
#define IN_VCF (1 << 0)  // 1
#define IN_DM (1 << 1)   // 2
#define IN_DXY (1 << 2)  // 4

#define A_BASE_IDX 0
#define C_BASE_IDX 1
#define G_BASE_IDX 2
#define T_BASE_IDX 3

/* ========================================================================== */

/// @brief strIsNumeric - check if string is numeric
/// @param val          - string to be checked
/// @return             - 1 if string is numeric, 0 otherwise
int strIsNumeric(const char *val);

// print generic usage information
void print_help(FILE *fp);

/// print formula usage information; to be used in formula specific errors
void print_help_formula(FILE *fp);

// verbosity level external global variable
// -v 0 or none -> 0000 0000 -> 0 -> verbose off [default, set in argStruct.cpp]
// -v 1 -> 0000 0001 -> 1 -> set the first bit, 0-indexed (1-1=0) verbose on with level 1
// -v 2 -> 0000 0010 -> 2 -> set the second bit, 0-indexed (2-1=1) verbose on with level 2
// ... and so on, up to -v 8
// -v 8 -> 1000 0000 -> 128 -> set the 8th bit, 0-indexed (8-1=7) verbose on with level 8
extern u_char VERBOSE;

const double NEG_INF = -std::numeric_limits<double>::infinity();

int extractDigits(int num, int digits);

char *get_time();

#endif
