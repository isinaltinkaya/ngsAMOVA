#ifndef __IO__
#define __IO__

#include <htslib/tbx.h>
#include <zlib.h>

#include "argStruct.h"

struct pairStruct;

typedef gzFile_s GZFILE;

typedef struct outfile_t outfile_t;
struct outfile_t {

    uint8_t ctype;    // output file compression type (e.g. none, gz, bgz. see PROGRAM_OUTFILE_CTYPE_*

    char* fn;        // full output filename (prefix + suffix + extension + cextension)
    kstring_t kbuf;  // kstring buffer for output file content

    void* fp;        // file pointer (FILE*, GZFILE*, BGZF*)
};

outfile_t* outfile_init(const char* suffix, const char* extension, const uint8_t ctype);
void outfile_destroy(outfile_t* outfile);
void outfile_write(outfile_t* outfile);
void outfile_write(outfile_t* outfile, const char* description);



//********************************************************************************
//*****************************  IO  *********************************************
//********************************************************************************

namespace IO {

    // /// @brief verbose - check if verbose level meets the threshold
    // ///
    // /// @param verbose_threshold - threshold to check against the verbose arg value
    // /// @return 1 if verbose level meets the threshold, 0 otherwise
    // ///
    // /// @example
    // /// if `-v 3` is used; (sets the bit at index 2 (3-1==2)
    // /// verbose(2) will check if a bit with index 1 (2-1==1) or higher is set
    // /// if(verbose(2)) will run;
    // ///
    // /// if `-v 1` is used; (sets the bit at index 0 (1-1==0)
    // /// verbose(2) will check if a bit with index 1 (2-1==1) or higher is set
    // /// if(verbose(2)) will not run;
    // int verbose(const int verbose_threshold);

    // TODO use a global count for count write warnings to log for x amount of sites like in angsd 
    // void vprint_count(const char* format, ...);

    /// @brief vprint - verbose print: print to stderr if verbose > 0
    /// @example vprint("Hello %s", "World");
    void vprint(const char* format, ...);

    /// @brief <o> verbose print with threshold
    /// @param verbose_threshold threshold for printing the specified message
    /// @example vprint(1, "Hello %s", "World"); // will print if verbose >= 1
    void vprint(const int verbose_threshold, const char* format, ...);

    /// @brief <o> verbose print with threshold and file pointer
    /// @param fp file pointer to print to
    /// @param verbose_threshold threshold for printing the specified message
    /// @example vprint(1, "Hello %s", "World"); // will print to fp if verbose >= 1
    void vprint(FILE* fp, const int verbose_threshold, const char* format, ...);

    hts_idx_t* load_bcf_csi_idx(const char* fn);
    tbx_t* load_vcf_tabix_idx(const char* fn);

    void requireArgFile(const char* fn, const char* requiredArg, const char* requiredFor);
    void requireArgFile(const char* fn, const char* requiredArg);

    void requireArgStr(const char* str, const char* requiredArg, const char* requiredFor);
    void requireArgStr(const char* str, const char* requiredArg);

    extern const char* FILE_EXTENSIONS[];

    int fileExists(const char* fn);

    char* setFileName(const char* a, const char* b);
    char* setFileName(const char* fn, const char* suffix, const char* fc_ext);

    const char* getFileExtension(const char* fn);

    FILE* getFile(const char* fname, const char* mode);
    gzFile getGzFile(const char* fname, const char* mode);
    FILE* openFileW(const char* a, const char* b);
    FILE* openFileW(char* c);
    gzFile openGzFileW(const char* a, const char* b);
    gzFile openGzFileW(char* c);

    bool isGzFile(const char* fn);

    namespace readFile {
        char* getFirstLine(FILE* fp);


    };  // namespace readFile

    namespace readGzFile {
        char* getFirstLine(char* fn);
        char* getFirstLine(gzFile fp);

    };  // namespace readGzFile

    namespace inspectFile {
        int count_nCols(const char* line, const char* delims);
        int count_nRows(char* fn, int HAS_COLNAMES);
        int count_nRows(FILE* fp, int HAS_COLNAMES);
    };  // namespace inspectFile

}  // namespace IO

kstring_t* kbuf_init();
void kbuf_destroy(kstring_t* kbuf);

#endif  // __IO__