#ifndef __IO__
#define __IO__

#include <zlib.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "argStruct.h"

struct sampleStruct;
struct formulaStruct;
struct pairStruct;


//********************************************************************************
//*****************************  IO  *********************************************
//********************************************************************************

namespace IO
{

    extern const char* FILE_EXTENSIONS[];

    char *setFileName(const char *a, const char *b);
    char *setFileName(const char *fn, const char *suffix, const char *fc_ext);

    const char *getFileExtension(const char *fn);

    FILE *getFile(const char *fname, const char *mode);
    gzFile getGzFile(const char *fname, const char *mode);
    FILE *openFileW(const char *a, const char *b);
    FILE *openFileW(char *c);
    gzFile openGzFileW(const char *a, const char *b);
    gzFile openGzFileW(char *c);

    int isGzFile(const char *fn);

    namespace readFile
    {
        int SFS(FILE *in_sfs_fp, const char *delims, sampleStruct *SAMPLES);

        char *getFirstLine(char *fn);
        char *getFirstLine(FILE *fp);

        int readToBuffer(char *fn, char **buffer_p, size_t *buf_size_p);

        int getBufferSize(FILE *fp);
        int getBufferSize(char *fn);
    };

    namespace readGzFile
    {
        char *getFirstLine(char *fn);
        char *getFirstLine(gzFile fp);

        int readToBuffer(char *fn, char **buffer_p, size_t *buf_size_p);
        // int getBufferSize(gzFile *fp);
        // int getBufferSize(char *fn);
    };

    namespace inspectFile
    {
        int count_nColumns(char *line, const char *delims);
        int count_nRows(char *fn, int HAS_COLNAMES);
        int count_nRows(FILE *fp, int HAS_COLNAMES);
    };

    namespace validateFile
    {
        int Metadata(FILE *in_mtd_fp, int nInds, int *keyCols, formulaStruct *FORMULA, const char *delims, int HAS_COLNAMES);
    };

    typedef struct outputStruct
    {

        char *fn = NULL;

        OUTFC fc;

        FILE *fp = NULL;
        gzFile gzfp = NULL;
        BGZF *bgzfp = NULL;
        // int print_header = 0; //TODO

        outputStruct(const char *fn_, const char *suffix, int fc_)
        {
            fc = OUTFC(fc_);
            switch (fc)
            {
            case OUTFC::NONE:
                fn = setFileName(fn_, suffix, FILE_EXTENSIONS[fc]);
                fp = openFileW(fn);
                break;
            case OUTFC::GZ:
                NEVER();
                fn = setFileName(fn_, suffix, FILE_EXTENSIONS[fc]);
                gzfp = openGzFileW(fn);
                break;
            case OUTFC::BBGZ:
                fn = setFileName(fn_, suffix, FILE_EXTENSIONS[fc]);
                bgzfp = bgzf_open(fn, "wb");
                break;
            default:
                fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
                exit(1);
                break;
            }
            fprintf(stderr,"\n[INFO] Opening output file: %s with compression type: %d (%s)\n", fn, fc, OUTFC_LUT[(OUTFC)fc]);

        }

        ~outputStruct()
        {
            // flush();
            
            fprintf(stderr,"\n[INFO] Closing output file: %s with compression type: %d (%s)\n", fn, fc, OUTFC_LUT[(OUTFC)fc]);

            switch (fc)
            {
            case OUTFC::NONE:
                FCLOSE(fp);
                break;
            case OUTFC::GZ:
                GZCLOSE(gzfp);
                break;
            case OUTFC::BBGZ:
                BGZCLOSE(bgzfp);
                break;
            default:
                fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
                exit(1);
                break;
            }
            FREE(fn);
        }

        void flush()
        {
            switch (fc)
            {
            case OUTFC::NONE:
                fflush(fp);
                break;
            case OUTFC::GZ:
                gzflush(gzfp, Z_SYNC_FLUSH);
                break;
            case OUTFC::BBGZ:
                ASSERT(bgzf_flush(bgzfp) == 0);
                break;
            default:
                fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
                exit(1);
                break;
            }
        }

        void *get_fp()
        {
            switch (fc)
            {
            case OUTFC::NONE:
                return fp;
                break;
            case OUTFC::GZ:
                return gzfp;
                break;
            case OUTFC::BBGZ:
                return bgzfp;
                break;
            default:
                fprintf(stderr, "\n[ERROR] Unknown file compression type (%d)\n", fc);
                exit(1);
                break;
            }
        }

        void write(const char *buf);
        void write(const kstring_t *kbuf);

    } outputStruct;

    typedef struct outFilesStruct
    {
        outputStruct *out_dm_fs = NULL;
        outputStruct *out_amova_fs = NULL;
        outputStruct *out_dev_fs = NULL;
        outputStruct *out_jgcd_fs = NULL;
        outputStruct *out_jgpd_fs = NULL;

        outFilesStruct(argStruct *args);

        ~outFilesStruct();

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
        // void Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, pairStruct *pair, argStruct *args, const char *sample1, const char *sample2);

        // void Sfs(const char *TYPE, IO::outputStruct *out_sfs_fs, argStruct *args, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);

        void M_PWD(const char *TYPE, IO::outputStruct *out_dm_fs, int nIndCmb, double *M_PWD);

    }

}
kstring_t *kbuf_init();
void kbuf_destroy(kstring_t *kbuf);

#endif // __IO__