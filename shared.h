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

#define ARG_DOEM_3GL 1
#define ARG_DOEM_10GL 2



#define PROGRAM_NEEDS_INDNAMES \
    ( ((0 != (args->doPhylo))))

#define PROGRAM_NEEDS_METADATA \
    ( ( args->doAMOVA || args->doDxy ))

#define PROGRAM_NEEDS_FORMULA \
    ( ( args->doAMOVA || args->doDxy ))

#define ARG_DOAMOVA_UNSET 0
#define ARG_DOAMOVA_SINGLERUN 1
#define ARG_DOAMOVA_BOOTSTRAP 2
#define ARG_DOAMOVA_PERMTEST 3


#define PROGRAM_WILL_PERFORM_BLOCK_BOOTSTRAPPING \
    ( (args->doAMOVA==ARG_DOAMOVA_BOOTSTRAP))
//TODO add block bootstrap test for nj and dxy



/// ----------------------------------------------------------------------- ///
// ARGUMENT VALUES

#define ARG_DOJGTM_UNSET 0
#define ARG_DOJGTM_3GT 1
#define ARG_DOJGTM_10GT 2

#define ARG_INTPLUS_INPUT_UNSET      (0<<0)
#define ARG_INTPLUS_INPUT_VCF        (1<<0)
#define ARG_INTPLUS_INPUT_DM         (1<<1)
#define ARG_INTPLUS_INPUT_MULTIDM    (1<<2)
#define ARG_INTPLUS_INPUT_DXY        (1<<3)
#define ARG_INTPLUS_INPUT_METADATA   (1<<4)
#define ARG_INTPLUS_INPUT_MAJORMINOR (1<<5)
#define ARG_INTPLUS_INPUT_ANCDER     (1<<6)
#define ARG_INTPLUS_INPUT_BLOCKS     (1<<7)
#define ARG_INTPLUS_INPUT_REGIONS    (1<<8)

#define PROGRAM_HAS_INPUT_VCF \
    ( (args->in_ft & ARG_INTPLUS_INPUT_VCF) )

#define PROGRAM_HAS_INPUT_DM \
    ( (args->in_ft & ARG_INTPLUS_INPUT_DM) )

#define PROGRAM_HAS_INPUT_MULTIDM \
    ( (args->in_ft & ARG_INTPLUS_INPUT_MULTIDM) )

#define PROGRAM_HAS_INPUT_DXY \
    ( (args->in_ft & ARG_INTPLUS_INPUT_DXY) )

#define PROGRAM_HAS_INPUT_METADATA \
    ( (args->in_ft & ARG_INTPLUS_INPUT_METADATA) )

#define PROGRAM_HAS_INPUT_MAJORMINOR \
    ( (args->in_ft & ARG_INTPLUS_INPUT_MAJORMINOR) )

#define PROGRAM_HAS_INPUT_ANCDER \
    ( (args->in_ft & ARG_INTPLUS_INPUT_ANCDER) )

#define PROGRAM_HAS_INPUT_BLOCKS \
    ( (args->in_ft & ARG_INTPLUS_INPUT_BLOCKS) )

#define PROGRAM_HAS_INPUT_REGIONS \
    ( (args->in_ft & ARG_INTPLUS_INPUT_REGIONS) )




#define ARG_DOMAJORMINOR_UNSET (0)
#define ARG_DOMAJORMINOR_BCF_REFALT1 (1)
#define ARG_DOMAJORMINOR_INFILE (2)

#define PROGRAM_WILL_USE_ALLELES_REF_ALT1 \
    ( ((args->doMajorMinor) == ARG_DOMAJORMINOR_BCF_REFALT1) )


/// ----------------------------------------------------------------------- ///
// BCF DATA SOURCE TAGS

#define ARG_INTPLUS_BCFSRC_FMT_GL (1<<0)
#define ARG_INTPLUS_BCFSRC_FMT_GT (1<<1)


#define PROGRAM_WILL_USE_BCF_FMT_GL \
    ( (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GL) )

#define PROGRAM_WILL_USE_BCF_FMT_GT \
    ( (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GT) )

#define PROGRAM_WILL_USE_RNG \
    ( (args->doAMOVA==ARG_DOAMOVA_BOOTSTRAP) || (args->doAMOVA==ARG_DOAMOVA_PERMTEST) )

/* ========================================================================== */
/* MACRO DEFINITIONS ======================================================== */
/* ========================================================================== */

/* -> FUNCTION-LIKE MACROS ---------------------------------------------------*/

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))



/// @brief MATRIX_GET_INDEX_* - get the index of the given i,j pair in the given * type matrix
///
/// LTED: assume i>j
///       for(i=1;i<n;++i) for(j=0;j<i;++j)
#define MATRIX_GET_INDEX_LTED_IJ(i,j) (((i)*((i)-1)/2+(j)))
#define MATRIX_GET_INDEX_LTED_IJ_UNORDERED(i,j) (((i)>(j)) ? MATRIX_GET_INDEX_LTED_IJ(i,j) : MATRIX_GET_INDEX_LTED_IJ(j,i))

/// LTID: assume i<=j
///       for(i=0;i<n;++i) for(j=0;j<=i;++j)
#define MATRIX_GET_INDEX_LTID_IJ(i,j) (((i)*((i)+1)/2+(j)))

/// UTED: assume i<j
///       for(i=0;i<n-1;++i) for(j=i+1;j<n;++j)
#define MATRIX_GET_INDEX_UTED_IJ(i,j,n) ((((i)*(n))-((i)*((i)+3))/2+(j)-1))

/// UTID: assume i>=j
///       for(i=0;i<n;++i) for(j=i;j<n;++j)
#define MATRIX_GET_INDEX_UTID_IJ(i,j,n) ((((i)*(n))-((i)*((i)+1))/2+(j)))

/// @brief MATRIX_GET_IJ_FROM_LTED_INDEX - get the i,j pair from the given linear index in a UTED matrix
#define MATRIX_GET_IJ_FROM_UTED_INDEX(idx,n,i,j) \
do{ \
    (i) = (n) - 2 - floor(sqrt(-8 * (lidx) + 4 * (n) * ((n) - 1) - 7) / 2.0 - 0.5); \
    (j) = (idx) + (i) + 1 - (n) * ((n) - 1) / 2 + ((n) - (i)) * (((n) - (i)) - 1) / 2; \
}while(0); 



#define CHECK_ARG_INTERVAL_INT(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %d'] Allowed range is [%d,%d]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)


#define CHECK_ARG_INTERVAL_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)

#define CHECK_ARG_INTERVAL_IE_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) >= (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)


#define CHECK_ARG_INTERVAL_II_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)


#define CHECK_ARG_INTERVAL_EI_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) <= (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)


#define CHECK_ARG_INTERVAL_EE_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) <= (minval)) || ((argval) >= (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)



#define CHECK_ARG_INTERVAL_01(argval, argstr) \
do { \
    if ( ((argval)!=0) && ((argval)!=1) ) { \
        \
            ERROR("Argument %s with value %d is out of range. Allowed values are 0 (for on/enable) and 1 (for off/disable)", (argstr), (argval)); \
    } \
} while(0)





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
#define ERROR(...)  \
    do{ \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n*******\n"); \
        exit(1); \
    } while(0)

  /*
   * Macro:[NEVER]
   * indicates that a point in the code should never be reached
   */
#define NEVER \
    do{ \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
        fprintf(stderr, "Control should never reach this point; please report this to the developers."); \
        fprintf(stderr, "\n*******\n"); \
        exit(1); \
    } while(0)


   /*
    * Macro:[ASSERT]
    * evaluate an expression, works the same way as the C-macro assert
    * except that DEBUG does not affect it (it is always active)
    * also prints the file and line info and exits the program
    * if the expression evaluates to false
    */
#define ASSERT_EXPAND(x) x
#define ASSERT(expr)                                                        \
    do {                                                                    \
        if (!(ASSERT_EXPAND(expr))){                                        \
            fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", \
                    __FILE__, __FUNCTION__, __LINE__, #expr);               \
            exit(1);                                                        \
        }                                                                   \
    } while(0)


    /*
    * Macro:[WARN]
    * print a custom warning message
    */
#define WARN(...)                                                            \
do {                                                                     \
	fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, __FUNCTION__, \
			__LINE__);                                                   \
	fprintf(stderr, __VA_ARGS__);                                        \
	fprintf(stderr, "\n");                                               \
} while(0)

    /*
    * Macro:[VWARN]
    * print a custom warning message if verbose is set
    */
#define VWARN(...)                                                 \
do {                                                           \
	if (0 != VERBOSITY_LEVEL) {                                  \
		fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, \
				__FUNCTION__, __LINE__);                       \
		fprintf(stderr, __VA_ARGS__);                          \
		fprintf(stderr, "\n");                                 \
	}                                                          \
} while(0)

#define VRUN(...) \
do { \
    if (0 != VERBOSITY_LEVEL) { \
        __VA_ARGS__; \
    } \
} while(0)




    //TODO also write to argsfile
         /*
          * Macro:[LOG]
          */
#define LOG(...) \
do{ \
    fprintf(stderr, "\n[LOG](%s)\t", __FUNCTION__); \
    fprintf(stderr, __VA_ARGS__); \
} while(0)

#define LOGADD(...) \
do{ \
    fprintf(stderr, "\n[LOG]\t"); \
    fprintf(stderr, __VA_ARGS__); \
} while(0)


          /*
           * Macro:[BEGIN_LOGSECTION]
           * print a custom message to stderr
           * indicating the start of a new section in the log
           */

#define BEGIN_LOGSECTION \
do{ \
    fprintf(stderr, "\n________________________________________\n"); \
    fprintf(stderr, "[LOG]\t-> Program is now at <%s/%s>\n", __FILE__, __FUNCTION__); \
}while(0)

#define END_LOGSECTION \
do{ \
    fprintf(stderr, "\n________________________________________\n"); \
}while(0)

#define BEGIN_LOGSECTION_MSG(msg) \
do{ \
    fprintf(stderr, "\n________________________________________\n"); \
    fprintf(stderr, "[LOG]\t-> Program is now at <%s/%s>\n", __FILE__, __FUNCTION__); \
    fprintf(stderr, "[LOG]\t-> %s\n", msg); \
}while(0)



           /*
           * Macro:[ASSERTM]
           * shortcut to evaluate an expression, works the same way as ASSERT macro
           * but also prints a custom message
           */
#define ASSERTM(expr, msg)                                                                               \
do{ \
    if (!(expr)) {                                                                                       \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n", __FILE__, __FUNCTION__, __LINE__, #expr); \
        fprintf(stderr, "\n%s!!\n", #msg);                                                               \
        exit(1);                                                                                         \
    } \
} while(0)


           /*
            * Macro:[REALLOC] - safe realloc with auto dimension multiplication
            * 						i.e. newsize * (sizeof *ptr)
            *
            * check if reallocation was successful
            * on success, set original pointer to new memory
            * otherwise print error message and exit
            *
            *
            * e.g. int *p=(int*)malloc( 2*sizeof(int) );
            *
            *      then
            *
            * 		int *tmp = realloc(p, 10*sizeof(int) );
            * 		if(tmp!=NULL){
            * 			p=tmp;
            * 		}else{
            * 			exit(1);
            * 		}
            *
            * 		== becomes ==
            *
            * 		REALLOC(int*,p,10);
            *
            *
            */
#define REALLOC(type, p, size) \
	do{ \
		void* tmp = (void*) realloc(((p)), (((size)) * sizeof( *((p)) ))); \
        if (NULL == tmp) { \
            fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
            fprintf(stderr, "Failed to reallocate memory\n");                                        \
            fprintf(stderr, "\n*******\n");                                                          \
            exit(1);                                                                                 \
        } \
        (p) = (type) tmp; \
        tmp = NULL; \
    } while(0)


#define MALLOC(type, p, size)  \
	do{ \
        ((p)) = NULL; \
		((p)) = ( ((type)) ) malloc ( ((size)) * sizeof( *((p))) ); \
		if(NULL == ((p))) { \
			fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
			fprintf(stderr, "Failed to allocate memory\n");                                        \
			fprintf(stderr, "\n*******\n");                                                          \
			exit(1);                                                                                 \
		} \
	}while(0);

            /*
             * Macro:[FREE]
             * force free
             *
             * shortcut to free memory and set pointer to NULL
             * and catch double free
             * throw error if pointer is NULL
             */
#define FREE(x)                                 \
	do{ \
		if (NULL != x) {                            \
			free(x);                                \
			x = NULL;                               \
		} else {                                    \
			ERROR("Trying to free a NULL pointer"); \
		} \
	}while(0);



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


#define PCLOSE(fp) \
do{ \
    if(0!=pclose((fp))){ \
        fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #fp); \
        exit(1); \
    } \
}while(0);

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

#define MAX_N_GROUPS_PER_LEVEL 10

#define MAX_N_BITS 64

// maximum level of verbosity allowed
#define MAX_PROGRAM_VERBOSITY_LEVEL 3

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
// #define BUF_NSITES 4096
#define BUF_NSITES 500000
#define BUF_NCONTIGS 256
#define BUF_NTOTSITES 4096

#define FREAD_BUF_SIZE 4096

// TODO check all references
#define FGETS_BUF_SIZE 1024

#define NSITES_BUF_INIT 4096

/* CONSTANTS -----------------------------------------------------------------*/

#define PROGRAM_PLOIDY 2

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

/* ========================================================================== */

/// @brief strIsNumeric - check if string is numeric
/// @param val          - string to be checked
/// @return             - 1 if string is numeric, 0 otherwise
int strIsNumeric(const char* val);

// print generic usage information
void print_help(FILE* fp);

/// print formula usage information; to be used in formula specific errors
void print_help_formula(FILE* fp);

const double NEG_INF = -std::numeric_limits<double>::infinity();

int extractDigits(int num, int digits);

char* get_time();

#endif
