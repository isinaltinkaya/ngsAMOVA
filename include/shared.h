/**
 * @file    shared.h
 * @brief   shared header file
 * @details contains shared stuff for the entire program
 */
#ifndef __SHARED_H__
#define __SHARED_H__

/* INCLUDES ----------------------------------------------------------------- */
#include <htslib/vcf.h>

#include "dev.h"
#include "version.h"
#include "build.h"
/* END-OF-INCLUDES ---------------------------------------------------------- */



/* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct argStruct argStruct;
typedef struct paramStruct paramStruct;
/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */

// ARG_INTPLUS_UNSET - INT+ argument unset value
#define ARG_INTPLUS_UNSET (0)


// -------------------------------
// -> ARG_PRINT_* - print options
#define ARG_PRINT_NOPRINT (0)

// --print-jgtm 0: do not print jgtm
// --print-jgtm 1: print jgtm
#define ARG_PRINT_JGTM (1)

// --print-dm 0: do not print dmat
// --print-dm 1: print orig dmat
// --print-dm 2: print verbose dmat
// --print-dm 3: print orig dmat and verbose dmat
#define ARG_INTPLUS_PRINT_DM_ORIG  (1<<0)
#define ARG_INTPLUS_PRINT_DM_VERBOSE  (1<<1)

// --print-pruned-dm 0: do not print pruned dmat
// --print-pruned-dm 1: print orig pruned dmat 
// --print-pruned-dm 2: print verbose pruned dmat
#define ARG_INTPLUS_PRINT_PRUNED_DM_ORIG (1<<0)
#define ARG_INTPLUS_PRINT_PRUNED_DM_VERBOSE (1<<1)

// --print-amova 0: do not print amova
// --print-amova 1: print amova as csv
// --print-amova 2: print amova as table
#define ARG_INTPLUS_PRINT_AMOVA_CSV   (1<<0)
#define ARG_INTPLUS_PRINT_AMOVA_TABLE (1<<1)

// --print-blocks 0: do not print blocks
// --print-blocks 1: print blocks
#define ARG_PRINT_BLOCKS (1)

// --print-bootstrap 0: do not print bootstrap samples
// --print-bootstrap 1: print bootstrap samples
#define ARG_PRINT_BOOTSTRAP (1)

// --print-tree 0: do not print tree
// --print-tree 1: print tree
#define ARG_PRINT_TREE (1)

// --print-dxy 0: do not print dxy
// --print-dxy 1: print dxy
#define ARG_PRINT_DXY (1)

// --print-ibd 0: do not print ibd
// --print-ibd 1: print ibd
// --print-ibd 2: print ibd persite scores
// --print-ibd 3: print ibd persite smoothed scores
#define ARG_INTPLUS_PRINT_IBD_SEGMENTS  (1<<0)
#define ARG_INTPLUS_PRINT_IBD_PERSITE_IBD_SCORES  (1<<1)
#define ARG_INTPLUS_PRINT_IBD_PERSITE_SMOOTHED_IBD_SCORES  (1<<2)

// -------------------------------


#define ARG_INTPLUS_BCFSRC_FMT_GL (1<<0)
#define ARG_INTPLUS_BCFSRC_FMT_GT (1<<1)

#define ARG_DOAMOVA_UNSET (0)
#define ARG_DOAMOVA_SINGLERUN (1)
#define ARG_DOAMOVA_BOOTSTRAP (2)
#define ARG_DOAMOVA_PERMTEST (3)

#define ARG_DOJGTM_UNSET (0)
#define ARG_DOJGTM_3GT (1)
#define ARG_DOJGTM_10GT (2)

#define ARG_DODIST_UNSET (0)
#define ARG_DODIST_CALCULATE (1)
#define ARG_DODIST_INFILE (2)
#define ARG_DODIST_WITH_BOOTSTRAP (3)

#define ARG_DOIBD_UNSET (0)
// TODO may be unnecessary if we already give --bcf-src
#define ARG_DOIBD_GL_METHOD (1)
#define ARG_DOIBD_GT_METHOD (2)

#define ARG_INTPLUS_INPUT_UNSET      (0<<0)
#define ARG_INTPLUS_INPUT_VCF        (1<<0)
#define ARG_INTPLUS_INPUT_DM         (1<<1)
#define ARG_INTPLUS_INPUT_MULTIDM    (1<<2)
#define ARG_INTPLUS_INPUT_METADATA   (1<<3)
#define ARG_INTPLUS_INPUT_MAJORMINOR (1<<4)
#define ARG_INTPLUS_INPUT_BLOCKS     (1<<6)
#define ARG_INTPLUS_INPUT_REGIONS    (1<<7)
#define ARG_INTPLUS_INPUT_JGTMAT     (1<<8)

#define ARG_DOMAJORMINOR_UNSET (0)
#define ARG_DOMAJORMINOR_BCF_REFALT1 (1)
#define ARG_DOMAJORMINOR_INFILE (2)
#define ARG_DOMAJORMINOR_BCF_AC (3)


// TODO check if thisis used currently
// 0 = basic (uninformed)
#define ARG_ALLOC_STRATEGY_BASIC (0)

// 1 = informed (init_nsites=vcfd->maxnsites/2, step_size=(vcfd->maxnsites)/8)
#define ARG_ALLOC_STRATEGY_INIT_MAXDIV2_STEP_MAXDIV8 1

// 2 = fully informed (init to vcfd->maxnsites, step size = 1024)
#define ARG_ALLOC_STRATEGY_INIT_MAX_STEP_1024 2


/*
 * Macro:[AT]
 * inject the file and line info as string
 */
#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)



 /// @brief PROGRAM_CURRENT_FUNCTION - get the current function name
 /// @note  use instead of __FUNCTION__

 /// source: https://www.boost.org/doc/libs/1_84_0/boost/current_function.hpp
#ifdef PROGRAM_DISABLE_CURRENT_FUNCTION
# define PROGRAM_CURRENT_FUNCTION " "

#elif defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) || (defined(__ICC) && (__ICC >= 600)) || defined(__ghs__) || defined(__clang__)
// # define PROGRAM_CURRENT_FUNCTION __PRETTY_FUNCTION__
# define PROGRAM_CURRENT_FUNCTION __FUNCTION__ // prefer non-pretty function

#elif defined(__DMC__) && (__DMC__ >= 0x810)
// # define PROGRAM_CURRENT_FUNCTION __PRETTY_FUNCTION__
# define PROGRAM_CURRENT_FUNCTION __FUNCTION__ // prefer non-pretty function

#elif defined(__FUNCSIG__)
# define PROGRAM_CURRENT_FUNCTION __FUNCSIG__

#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))
# define PROGRAM_CURRENT_FUNCTION __FUNCTION__

#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)
# define PROGRAM_CURRENT_FUNCTION __FUNC__

#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
# define PROGRAM_CURRENT_FUNCTION __func__

#elif defined(__cplusplus) && (__cplusplus >= 201103)
# define PROGRAM_CURRENT_FUNCTION __func__

#else
# define PROGRAM_CURRENT_FUNCTION " "

#endif





 /*
  * Macro:[ERROR]
  * print a custom error message and exit the program
  */
#define ERROR(...)  \
    do{ \
        fprintf(stderr, "\n\n*******\n[ERROR](%s:%d)\n\t", __FILE__, __LINE__); \
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
        ERROR("Control should never reach this point; please report this to the developers."); \
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
                    __FILE__, PROGRAM_CURRENT_FUNCTION, __LINE__, #expr);               \
            exit(1);                                                        \
        }                                                                   \
    } while(0)


    /*
    * Macro:[WARN]
    * print a custom warning message
    */
#define WARN(...)                                                            \
do {                                                                     \
	fprintf(stderr, "\n\n[WARNING](%s): ", PROGRAM_CURRENT_FUNCTION); \
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
		fprintf(stderr, "\n\n[WARNING](%s): ", PROGRAM_CURRENT_FUNCTION); \
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
    fprintf(stderr, "\n[LOG]\t"); \
    fprintf(stderr, __VA_ARGS__); \
} while(0)

#define LONGLOG(...) \
do{ \
    fprintf(stderr, "\n[LOG](%s)\t", PROGRAM_CURRENT_FUNCTION); \
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
    fprintf(stderr, "[LOG]\t-> Program is now at: <%s>\n", PROGRAM_CURRENT_FUNCTION); \
}while(0)


#define BEGIN_LOGSECTION_MSG(...) \
do{ \
    fprintf(stderr, "\n________________________________________\n"); \
    fprintf(stderr, "[LOG]\t-> Program is now at: "); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, " <%s>\n", PROGRAM_CURRENT_FUNCTION); \
}while(0)

#define END_LOGSECTION \
do{ \
    fprintf(stderr, "\n\n[LOG]\t-> Finished: <%s>\n", PROGRAM_CURRENT_FUNCTION); \
    fprintf(stderr, "________________________________________\n"); \
}while(0)


#define END_LOGSECTION_MSG(...) \
do{ \
    fprintf(stderr, "\n\n[LOG]\t-> Finished: "); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, " <%s>\n", PROGRAM_CURRENT_FUNCTION); \
    fprintf(stderr, "________________________________________\n"); \
}while(0)


           /*
           * Macro:[ASSERTM]
           * shortcut to evaluate an expression, works the same way as ASSERT macro
           * but also prints a custom message
           */
#define ASSERTM(expr, msg)                                                                               \
do{ \
    if (!(expr)) {                                                                                       \
        fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n", __FILE__, PROGRAM_CURRENT_FUNCTION, __LINE__, #expr); \
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
            fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, PROGRAM_CURRENT_FUNCTION, __LINE__); \
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
			fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, PROGRAM_CURRENT_FUNCTION, __LINE__); \
			fprintf(stderr, "Failed to allocate memory\n");                                        \
			fprintf(stderr, "\n*******\n");                                                          \
			exit(1);                                                                                 \
		} \
	}while(0)

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
	}while(0)



             /*
             * Macro:[FCLOSE]
             * safe fclose
             */
#define FCLOSE(fp)             \
    do{ \
        if(NULL==(fp)){        \
            ERROR("Trying to fclose a NULL file pointer"); \
        } \
        if (fclose((FILE*)(fp)) != 0) { \
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #fp); \
            exit(1);           \
        }                      \
        (fp) = NULL;           \
    } while(0)


             /*
             * Macro:[BGZCLOSE]
             * safe bgzf_close
             */
#define BGZCLOSE(fp) \
do{ \
    if (NULL != (fp)) { \
        if (bgzf_close((BGZF*)(fp)) != 0) { \
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #fp); \
            exit(1); \
        } \
        (fp) = NULL; \
    } else { \
        ERROR("Trying to bgzf_close a NULL file pointer"); \
    } \
}while(0)

             /*
             * Macro:[GZCLOSE]
             * safe gzclose
             */
#define GZCLOSE(fp) \
do{ \
    if (NULL != (fp)) { \
        if (gzclose((gzFile_s*) (fp)) != Z_OK) { \
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #fp); \
            exit(1); \
        } \
        (fp) = NULL; \
    } else { \
        ERROR("Trying to gzclose a NULL file pointer"); \
    } \
}while(0)


             /*
             * Macro:[PCLOSE]
             * safe pclose
             */
#define PCLOSE(fp) \
do{ \
    if (NULL != (fp)) { \
        if (pclose((FILE*)(fp)) != 0) { \
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", __FILE__, __LINE__, #fp); \
            exit(1); \
        } \
        (fp) = NULL; \
    } else { \
        ERROR("Trying to pclose a NULL file pointer"); \
    } \
}while(0)

             /* LIMIT DEFINING MACROS -----------------------------------------------------*/

             // maximum number of tokens allowed in amova formula string
#define MAX_N_FORMULA_TOKENS 10

#define MAX_N_GROUPS_PER_LEVEL 10

#define MAX_N_BITS 64

// maximum level of verbosity allowed
#define MAX_PROGRAM_VERBOSITY_LEVEL 3


#define PROGRAM_OUTFILE_CTYPE_NONE (0)
#define PROGRAM_OUTFILE_CTYPE_RAW (1)
#define PROGRAM_OUTFILE_CTYPE_GZ (2)
#define PROGRAM_OUTFILE_CTYPE_BGZ (3)


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

// TODO check all references
#define FGETS_BUF_SIZE 1024

// -> PROGRAM_NEEDS_* - flags for checking if a specific things is needed for the program to run
// useful to check if the user has provided all necessary inputs before running the program

#define PROGRAM_NEEDS_INDNAMES \
    ( ((0 != (args->doPhylo))))

#define PROGRAM_NEEDS_METADATA \
    ( ( args->doAMOVA || args->doDxy ))

#define PROGRAM_NEEDS_FORMULA \
    ( ( args->doAMOVA || args->doDxy ))

#define PROGRAM_RUN_IS_DRYRUN \
    ( (args->doDryRun != 0) )

/// whether program needs full gldata to be read into memory or it can stream the sites (either one by one or in buffers)
#define PROGRAM_NEEDS_FULL_GLDATA \
    ( (PROGRAM_WILL_USE_BCF_FMT_GL) && (args->doEM) )

// -> PROGRAM_WILL - flags for checking if the program will perform a specific action

#define PROGRAM_WILL_USE_ALLELES_REF_ALT1 \
    ( ((args->doMajorMinor) == ARG_DOMAJORMINOR_BCF_REFALT1) )


//TODO add block bootstrap test for nj and dxy

#define PROGRAM_WILL_USE_BCF_FMT_GL \
    ( (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GL) )

#define PROGRAM_WILL_USE_BCF_FMT_GT \
    ( (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GT) )

#define PROGRAM_WILL_USE_RNG \
    ( (args->doAMOVA==ARG_DOAMOVA_BOOTSTRAP) || (args->doAMOVA==ARG_DOAMOVA_PERMTEST) )

// -> PROGRAM_HAS - flags for checking if the program has a specific input

#define PROGRAM_HAS_INPUT_VCF \
    ( (args->in_ft & ARG_INTPLUS_INPUT_VCF) )

#define PROGRAM_HAS_INPUT_DM \
    ( (args->in_ft & ARG_INTPLUS_INPUT_DM) )

#define PROGRAM_HAS_INPUT_JGTM \
    ( (args->in_ft & ARG_INTPLUS_INPUT_JGTMAT) )

#define PROGRAM_HAS_INPUT_MULTIDM \
    ( (args->in_ft & ARG_INTPLUS_INPUT_MULTIDM) )

#define PROGRAM_HAS_INPUT_DXY \
    ( (args->in_ft & ARG_INTPLUS_INPUT_DXY) )

#define PROGRAM_HAS_INPUT_METADATA \
    ( (args->in_ft & ARG_INTPLUS_INPUT_METADATA) )

#define PROGRAM_HAS_INPUT_MAJORMINOR \
    ( (args->in_ft & ARG_INTPLUS_INPUT_MAJORMINOR) )

#define PROGRAM_HAS_INPUT_BLOCKS \
    ( (args->in_ft & ARG_INTPLUS_INPUT_BLOCKS) )

#define PROGRAM_HAS_INPUT_REGIONS \
    ( (args->in_ft & ARG_INTPLUS_INPUT_REGIONS) )



// -> CHECK_ARG_* - argument checking macros

#define CHECK_ARG_INTERVAL_INT(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %d'] Allowed range is [%d,%d]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)


#define CHECK_ARG_INTERVAL_II_INT(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %d'] Allowed range is [%d,%d]", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)

#define CHECK_ARG_INTERVAL_EE_INT(argval, minval, maxval, argstr) \
do { \
    if (((argval) <= (minval)) || ((argval) >= (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %d'] Allowed range is (%d,%d)", (argstr), (argval), (minval), (maxval)); \
    } \
} while(0)


#define CHECK_ARG_INTERVAL_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
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


// for INT+ arguments, check if the value is within the range [0, (2^0+2^1+2^2+...+2^maxval]
#define CHECK_ARG_INTERVAL_INTPLUS(argval, maxval, argstr) \
do { \
    if ( ((argval) < 0) || ((argval) > ((1 << (maxval+1)) - 1)) ) { \
        \
            ERROR("Argument %s with value %d is out of range. Allowed values are in the range [0, %d]", (argstr), (argval), ((1 << (maxval+1)) - 1)); \
    } \
} while(0)






// -> other function-like macros


#define BYTES_TO_MB_DBL(x) ((double) ((x) / (1048576)))
#define BYTES_TO_GB_DBL(x) ((double) ((x) / (1073741824)))

#define GET_ARR_SIZE_BYTES_1D(type, n1) ( ((sizeof(type))*(n1) + 8))

#define GET_ARR_SIZE_BYTES_2D(type, n1, n2) ( (n1) * ( ((sizeof(type))*(n2) + 8)))

#define GET_ARR_SIZE_BYTES_3D(type, n1, n2, n3) ( (n1) * ( (sizeof(type)) * (n2) * (n3) + (8 * (n2)) + 8))

#define KS_INIT { 0, 0, NULL }

#define SQUARE(n) ((((n)) * ((n))))

#define POW10(x) (exp((x) * M_LN10))


/*
 * Macro:[LOG2LN]
 * convert log_10(x) to log_e(x)
 */
#define LOG2LN(x) ((((x)) / (M_LOG10E)))

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define ABS(x) (((x) < 0) ? -(x) : (x))

 /// @brief MATRIX_GET_INDEX_* - get the index of the given i,j pair in the given * type matrix
 ///
 /// LTED: assume i>j
 ///       for(i=1;i<n;++i) for(j=0;j<i;++j)
#define MATRIX_GET_INDEX_LTED_IJ(i,j) (((i)*((i)-1)/2+(j)))
#define MATRIX_GET_INDEX_LTED_IJ_UNORDERED(i,j) (((i)>(j)) ? MATRIX_GET_INDEX_LTED_IJ((i),(j)) : MATRIX_GET_INDEX_LTED_IJ((j),(i)))

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
}while(0) 


/* CONSTANTS -----------------------------------------------------------------*/


#define PROGRAM_PLOIDY 2

// 1/9
#define FRAC_1_9 0.1111111111111111


/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEFS ---------------------------------------------------------- */
/* END-OF-TYPEDEFS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */
/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */



/* ========================================================================== */
/* ENUMERATIONS ============================================================= */

// output file compression types
enum OUTFC {
    NONE,  // no compression [0]
    GZ,    // gzip [1]
    BBGZ,  // bgzip (binary) [2]
};

/* ========================================================================== */

#endif // __SHARED_H__