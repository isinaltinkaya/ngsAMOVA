#ifndef __IO__
#define __IO__

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <zlib.h>

#include "argStruct.h"

struct pairStruct;

//********************************************************************************
//*****************************  IO  *********************************************
//********************************************************************************

namespace IO {

void append_argsFile(const char *format, ...);

/// @brief validateString - check if string is valid
void validateString(const char *str);

/// @brief verbose - check if verbose level meets the threshold
///
/// @param verbose_threshold - threshold to check against the verbose arg value
/// @return 1 if verbose level meets the threshold, 0 otherwise
///
/// @example
/// if `-v 3` is used; (sets the bit at index 2 (3-1==2)
/// verbose(2) will check if a bit with index 1 (2-1==1) or higher is set
/// if(verbose(2)) will run;
///
/// if `-v 1` is used; (sets the bit at index 0 (1-1==0)
/// verbose(2) will check if a bit with index 1 (2-1==1) or higher is set
/// if(verbose(2)) will not run;
int verbose(const int verbose_threshold);

/// @brief vprint - verbose print: print to stderr if verbose > 0
/// @example vprint("Hello %s", "World");
void vprint(const char *format, ...);

/// @brief <o> verbose print with threshold
/// @param verbose_threshold threshold for printing the specified message
/// @example vprint(1, "Hello %s", "World"); // will print if verbose >= 1
void vprint(const int verbose_threshold, const char *format, ...);

/// @brief <o> verbose print with threshold and file pointer
/// @param fp file pointer to print to
/// @param verbose_threshold threshold for printing the specified message
/// @example vprint(1, "Hello %s", "World"); // will print to fp if verbose >= 1
void vprint(FILE *fp, const int verbose_threshold, const char *format, ...);

/// @brief <o> verbose print with threshold and file pointer, prints to both stderr and file
/// @param fp file pointer to print to
/// @param verbose_threshold threshold for printing the specified message
/// @example vprint(1, "Hello %s", "World"); // will print to both stderr and fp if verbose >= 1
void vvprint(FILE *fp, const int verbose_threshold, const char *format, ...);

hts_idx_t *load_bcf_csi_idx(const char *fn);
tbx_t *load_vcf_tabix_idx(const char *fn);

void requireArgFile(const char *fn, const char *requiredArg, const char *requiredFor);
void requireArgFile(const char *fn, const char *requiredArg);

void requireArgStr(const char *str, const char *requiredArg, const char *requiredFor);
void requireArgStr(const char *str, const char *requiredArg);

extern const char *FILE_EXTENSIONS[];

int fileExists(const char *fn);

char *setFileName(const char *a, const char *b);
char *setFileName(const char *fn, const char *suffix, const char *fc_ext);

const char *getFileExtension(const char *fn);

FILE *getFile(const char *fname, const char *mode);
gzFile getGzFile(const char *fname, const char *mode);
FILE *openFileW(const char *a, const char *b);
FILE *openFileW(char *c);
gzFile openGzFileW(const char *a, const char *b);
gzFile openGzFileW(char *c);

bool isGzFile(const char *fn);

namespace readFile {
char *getFirstLine(const char *fn);
char *getFirstLine(FILE *fp);

char *readToBuffer(const char *fn);

};  // namespace readFile

namespace readGzFile {
char *getFirstLine(char *fn);
char *getFirstLine(gzFile fp);

int readToBuffer(char *fn, char **buffer_p, size_t *buf_size_p);
};  // namespace readGzFile

namespace inspectFile {
int count_nCols(const char *line, const char *delims);
int count_nRows(char *fn, int HAS_COLNAMES);
int count_nRows(FILE *fp, int HAS_COLNAMES);
};  // namespace inspectFile

/// @brief outputStruct - struct for output files
///
/// @param kbuf - kstring buffer for output file contents
///                 -> if kbuf will only be used once,
///                 memory is allocated and deallocated inplace
///                 (inside the specific writing function)
///                 e.g. see amovaStruct::print_as_csv()
///                 -> if kbuf will be used multiple times
///                 i.e. appending to kbuf from different functions
///                 memory is allocated and deallocated in the constructor
///                 and destructor respectively
///                 e.g. see out_args_fs
///                 for thread safety, use a separate kstring_t within each thread
///                 and append to the main kbuf after joining the threads
///
typedef struct outputStruct {
    char *fn = NULL;

    OUTFC fc;

    FILE *fp = NULL;
    gzFile gzfp = NULL;
    BGZF *bgzfp = NULL;

    kstring_t *kbuf = NULL;

    outputStruct(const char *fn_, const char *suffix, int fc_);
    ~outputStruct();

    void flush();

    void *get_fp();

    void write(const char *buf);
    void write(kstring_t *kbuf);
    void kbuf_write();  // write internal kbuf
    void kbuf_destroy_buffer();

} outputStruct;

typedef struct outFilesStruct {
    outputStruct *out_args_fs = NULL;
    outputStruct *out_dm_fs = NULL;
    outputStruct *out_amova_fs = NULL;
    outputStruct *out_dev_fs = NULL;
    outputStruct *out_jgcd_fs = NULL;
    outputStruct *out_dxy_fs = NULL;
    outputStruct *out_nj_fs = NULL;
    outputStruct *out_blockstab_fs = NULL;
    outputStruct *out_v_bootstrapRep_fs = NULL;

} outFilesStruct;

void outFilesStruct_init(outFilesStruct *ofs);
void outFilesStruct_destroy(outFilesStruct *ofs);

namespace print {

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

}  // namespace print

}  // namespace IO
extern IO::outFilesStruct *outFiles;

kstring_t *kbuf_init();
void kbuf_destroy(kstring_t *kbuf);

#endif  // __IO__